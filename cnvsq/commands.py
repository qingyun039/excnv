"""Command-line interface and corresponding API for CNVkit."""
# NB: argparse CLI definitions and API functions are interwoven:
#   "_cmd_*" handles I/O and arguments processing for the command
#   "do_*" runs the command's functionality as an API
import argparse
import logging
import os
import sys

import pandas as pd
from skgenome import tabio, GenomicArray
from cnvlib import core, segmentation, call
from cnvlib.cmdutil import (verify_sample_sex, read_cna)

from ._version import __version__
from . import reference, fix, params
from . import utils
from . import plotcnv
from .classifycnv import classifycnv


__all__ = []
def public(fn):
    __all__.append(fn.__name__)
    return fn


AP = argparse.ArgumentParser(
        description="CNVkit, a command-line toolkit for copy number analysis.",
        epilog="See the online manual for details: https://cnvkit.readthedocs.io")
AP_subparsers = AP.add_subparsers(
        required=True,
        help="Sub-commands (use with -h for more info)")


# _____________________________________________________________________________
# Core pipeline

# reference -------------------------------------------------------------------

do_reference = public(reference.do_reference)
do_reference_flat = public(reference.do_reference_flat)


def _cmd_reference(args):
    """Compile a coverage reference from the given files (normal samples)."""
    usage_err_msg = ("Give .bam samples OR .fa genome")
    mod_params(args)
    if args.bamfiles:
        # Pooled reference
        female_samples = ((args.sample_sex.lower() not in ['y', 'm', 'male'])
                          if args.sample_sex else None)
        ref_probes = reference.do_reference(args.bamfiles, args.binsize, args.fasta,
                                            args.male_reference, female_samples,
                                            args.do_gc, args.do_edge,
                                            args.do_rmask, args.cluster,
                                            args.min_cluster_size, args.bigwig,
                                            by_count=args.count,
                                            min_mapq=args.min_mapq,
                                            processes=args.processes)
    elif args.binsize:
        # Flat refence
        assert not args.fasta, usage_err_msg
        ref_probes = reference.do_reference_flat(args.binsize,
                                                 args.fasta,
                                                 args.bigwig,
                                                 args.male_reference)
    else:
        raise ValueError(usage_err_msg)

    ref_fname = args.output or "cnv_reference.cnn"
    core.ensure_path(ref_fname)
    tabio.write(ref_probes, ref_fname)


P_reference = AP_subparsers.add_parser('reference', help=_cmd_reference.__doc__)
P_reference.add_argument('bamfiles', nargs='*',
        help="""Normal-sample align .bam files, or the
                directory that contains them.""")
P_reference.add_argument('-b', '--binsize', type=int,
        help="Binsize, genome will splited (e.g. 10000)")
P_reference.add_argument('-f', '--fasta',
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)")
P_reference.add_argument('-B', '--bigwig',
        help="Alignability file, BigWig format (e.g. UCSC hg19.bw)")
P_reference.add_argument('-o', '--output', metavar="FILENAME",
        help="Output file name.")
P_reference.add_argument('-C', '--cluster',
        action='store_true',
        help="""Calculate and store summary stats for clustered subsets of the
                normal samples with similar coverage profiles.""")
P_reference.add_argument('--min-cluster-size',
        metavar="NUM",
        type=int,
        default=4,
        help="""Minimum cluster size to keep in reference profiles.""")
P_reference.add_argument('-x', '--sample-sex', '-g', '--gender',
        dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the chromosomal sex of all given samples as male or
                female. (Default: guess each sample from coverage of X and Y
                chromosomes).""")
P_reference.add_argument('-y', '--male-reference', '--haploid-x-reference',
        action='store_true',
        help="""Create a male reference: shift female samples' chrX
                log-coverage by -1, so the reference chrX average is -1.
                Otherwise, shift male samples' chrX by +1, so the reference chrX
                average is 0.""")

P_reference_bias = P_reference.add_argument_group(
    "To disable specific automatic bias corrections")
P_reference_bias.add_argument('--no-gc', dest='do_gc', action='store_false',
        help="Skip GC correction.")
P_reference_bias.add_argument('--do-edge', dest='do_edge', action='store_true',
        help="Skip edge-effect correction.")
P_reference_bias.add_argument('--no-rmask', dest='do_rmask', action='store_false',
        help="Skip RepeatMasker correction.")
P_reference.set_defaults(func=_cmd_reference)

P_coverage = P_reference.add_argument_group(
        "To load Coverage info from BAM file")
P_coverage.add_argument('-c', '--count', action='store_true',
        help="""Get read depths by counting read midpoints within each bin.
                (An alternative algorithm).""")
P_coverage.add_argument('-q', '--min-mapq', type=int, default=30,
        help="""Minimum mapping quality score (phred scale 0-60) to count a read
                for coverage depth.  [Default: %(default)s]""")
P_coverage.add_argument('-p', '--processes',
        nargs='?', type=int, const=0, default=1,
        help="""Number of subprocesses to calculate coverage in parallel.
                Without an argument, use the maximum number of available CPUs.
                [Default: use 1 process]""")

P_filter = P_reference.add_argument_group("Filter Bins")
P_filter.add_argument('--min-gc', type=float, help="Mininum GC content")
P_filter.add_argument('--max-gc', type=float, help="Maxinum GC content")
P_filter.add_argument('--min-alignability', type=float, help="Mininum Alignability")
P_filter.add_argument('--max-n-content', type=float, help="Maxinum N content")
P_filter.add_argument('--max-rmask', type=float, help="Maxinum RepeatMask content")
P_filter.add_argument('--max-zero-depth', type=float, help="Maxinum Zero depth fraction")
P_filter.add_argument('--min-ref-coverage', type=float, help="Mininum reference coverage")
P_filter.add_argument('--max-ref-spread', type=float, help="Mininum reference spread")


# fix -------------------------------------------------------------------------

do_fix = public(fix.do_fix)


def _cmd_detect(args):
    """Combine target and antitarget coverages and correct for biases.

    Adjust raw coverage data according to the given reference, correct potential
    biases and re-center.
    """
    mod_params(args)
    cnarr = fix.do_fix(args.bamfile, read_cna(args.reference),
                              args.do_gc, args.do_edge, args.do_rmask,
                              args.cluster,
                              by_count=args.count,
                              min_mapq=args.min_mapq,
                              processes=args.processes)
    if (not args.threshold) and (args.method == 'cbs'):
        args.threshold = 0.001
    cns = segmentation.do_segmentation(cnarr, args.method, args.threshold,
                                           skip_low=args.drop_low_coverage,
                                           skip_outliers=args.drop_outliers,
                                           rscript_path=args.rscript_path,
                                           processes=args.processes,
                                           smooth_cbs=args.smooth_cbs)
    if args.center_at:
        logging.info("Shifting log2 ratios by %f", -args.center_at)
        cns['log2'] -= args.center_at
    elif args.center:
        cns.center_all(args.center, skip_low=args.drop_low_coverage, verbose=True)

    is_sample_female = verify_sample_sex(cnarr, args.sample_sex, args.male_reference)
    result = call.do_call(cns, None, args.call_method, 2, 1, args.male_reference, is_sample_female, args.filters, args.thresholds)


    if args.output:
        outdir = os.path.dirname(args.output)
        if not outdir:
            outdir = './'
    else:
        outdir = './'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    tabio.write(cnarr, os.path.join(outdir, cnarr.sample_id + '.cnr'))
    tabio.write(result, os.path.join(outdir, cnarr.sample_id + '.cns'))
    
    result = utils.filter_normal(result, is_sample_female)
    result = utils.add_cnv_type(result, is_sample_female)
    tabio.write(result, args.output or os.path.join(outdir, cnarr.sample_id + '.tsv'))

    if args.plotcnv or args.annotation:
        otherdir = os.path.join(outdir, cnarr.sample_id)
        if not os.path.exists(otherdir):
            os.makedirs(otherdir)

    if args.plotcnv:
        if not args.ideogram:
            raise ValueError("plot CNV must specify idsogram file")
        ideo = plotcnv.load_ideo_data(args.ideogram)
        figout = os.path.join(otherdir, 'genome')
        plotcnv.genome_whole_view(ideo, result.data, min_cnv=10000, figout=figout, title=cnarr.sample_id)

        # chrom
        for chrom in cnarr.chromosome.unique():
            chrom_ideo = ideo[ideo.chrom == chrom]
            chrom_cnr = cnarr[cnarr.chromosome == chrom]
            chrom_cns = cns[cns.chromosome == chrom]
            # where female reference
            if chrom == 'chrY':
                if is_sample_female:
                    ref_copy = 0.1
                else:
                    ref_copy = 1
            else:
                ref_copy = 2
            figout = os.path.join(outdir, cnarr.sample_id, f'{chrom}')
            plotcnv.chromosome_cnv_view(chrom_ideo, chrom_cnr.data, chrom_cns.data, ref_copy=ref_copy, figout=figout, title=cnarr.sample_id + '_' + chrom)

    if args.annotation:
        cnv_list = set(result.data[['chromosome', 'start', 'end', 'type']]
               .agg(lambda x: '_'.join(x.values.astype(str)), axis=1)
               .T)
        outfile = os.path.join(otherdir + 'cna')
        args.GenomeBuild = 'hg19'
        args.precise = True
        args.outdir = otherdir
        args.cores = args.processes
        classifycnv.args = args
        classifycnv.do_classifycnv(cnv_list, outfile)




def csvstring(text):
    return tuple(map(float, text.split(",")))

P_detect = AP_subparsers.add_parser('detect', help=_cmd_detect.__doc__)
P_detect.add_argument('bamfile',
        help="sample align file (.bam).")
P_detect.add_argument('reference',
        help="Reference coverage (.cnn).")
P_detect.add_argument('-C', '--cluster',
        action='store_true',
        help="""Compare and use cluster-specific values present in the
                reference profile. (Requires that the reference profile
                was built with the --cluster option.)""")
P_detect.add_argument('-i', '--sample-id',
        help="Sample ID for target/antitarget files. Otherwise inferred from file names.")
# P_detect.add_argument('--do-gc', action='store_true', default=True,
#         help="Do GC correction.")
# P_detect.add_argument('--do-edge', action='store_true',
#         help="Do edge-effect correction.")
# P_detect.add_argument('--do-size', action='store_true',
#         help="Do interval-size correction.")
P_detect.add_argument('--no-gc', dest='do_gc', action='store_false',
        help="Skip GC correction.")
P_detect.add_argument('--do-edge', dest='do_edge', action='store_true',
        help="Skip edge-effect correction.")
P_detect.add_argument('--no-rmask', dest='do_rmask', action='store_false',
        help="Skip RepeatMasker correction.")
P_detect.add_argument('-o', '--output', metavar="FILENAME",
        help="Output file name.")
P_detect.add_argument('-p', '--processes',
        nargs='?', type=int, const=0, default=1,
        help="""Number of subprocesses to calculate coverage in parallel.
                Without an argument, use the maximum number of available CPUs.
                [Default: use 1 process]""")

P_coverage = P_detect.add_argument_group(
        "To load Coverage info from BAM file")
P_coverage.add_argument('-c', '--count', action='store_true',
        help="""Get read depths by counting read midpoints within each bin.
                (An alternative algorithm).""")
P_coverage.add_argument('-q', '--min-mapq', type=int, default=30,
        help="""Minimum mapping quality score (phred scale 0-60) to count a read
                for coverage depth.  [Default: %(default)s]""")


P_segment = P_detect.add_argument_group(
        "Segmentation arguments")
P_segment.add_argument('-m', '--method',
        choices=segmentation.SEGMENT_METHODS,
        default='cbs',
        help="""Segmentation method (see docs), or 'none' for chromosome
                arm-level averages as segments. [Default: %(default)s]""")
P_segment.add_argument('-t', '--threshold', type=float,
        help="""Significance threshold (p-value or FDR, depending on method) to
                accept breakpoints during segmentation.
                For HMM methods, this is the smoothing window size.""")
P_segment.add_argument("--drop-low-coverage", action='store_true',
        help="""Drop very-low-coverage bins before segmentation to avoid
                false-positive deletions in poor-quality tumor samples.""")
P_segment.add_argument("--drop-outliers", metavar="FACTOR",
        type=float, default=10,
        help="""Drop outlier bins more than this many multiples of the 95th
                quantile away from the average within a rolling window.
                Set to 0 for no outlier filtering.
                [Default: %(default)g]""")
P_segment.add_argument("--rscript-path", metavar="PATH", default="Rscript",
        help="""Path to the Rscript excecutable to use for running R code.
                Use this option to specify a non-default R installation.
                [Default: %(default)s]""")
P_segment.add_argument('--smooth-cbs', action='store_false',
        help="""Perform an additional smoothing before CBS segmentation, 
								which in some cases may increase the sensitivity. 
                                Used only for CBS method.[Default: true]""")
P_call = P_detect.add_argument_group("Call arguments")
P_call.add_argument('-M', '--call-method',
        choices=('threshold', 'clonal', 'none'), default='clonal',
        help="""Calling method. [Default: %(default)s]""")
P_call.add_argument("--center", nargs='?', const='median',
        choices=('mean', 'median', 'mode', 'biweight'),
        help="""Re-center the log2 ratio values using this estimator of the
                center or average value. ('median' if no argument given.)""")
P_call.add_argument("--center-at", type=float,
        help="""Subtract a constant number from all log2 ratios. For "manual"
                re-centering, in case the --center option gives unsatisfactory
                results.)""")
P_call.add_argument('--filter', action='append', default=[], dest='filters',
        choices=('ampdel', 'cn', 'ci', 'sem', # 'bic'
                ),
        help="""Merge segments flagged by the specified filter(s) with the
                adjacent segment(s).""")
P_call.add_argument('-T', '--thresholds',
        type=csvstring, default="-1.1,-0.25,0.2,0.7",
        help="""Hard thresholds for calling each integer copy number, separated
                by commas. Use the '=' sign on the command line, e.g.: -t=-1,0,1
                [Default: %(default)s]""")
P_call.add_argument('-x', '--sample-sex', '-g', '--gender', dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the sample's chromosomal sex as male or female.
                (Otherwise guessed from X and Y coverage).""")
P_call.add_argument('-y', '--male-reference', '--haploid-x-reference',
        action='store_true',
        help="""Was a male reference used?  If so, expect half ploidy on
                chrX and chrY; otherwise, only chrY has half ploidy.  In CNVkit,
                if a male reference was used, the "neutral" copy number (ploidy)
                of chrX is 1; chrY is haploid for either reference sex.""")
P_filter = P_detect.add_argument_group("Filter Bins")
P_filter.add_argument('--min-gc', type=float, help="Mininum GC content")
P_filter.add_argument('--max-gc', type=float, help="Maxinum GC content")
P_filter.add_argument('--min-alignability', type=float, help="Mininum Alignability")
P_filter.add_argument('--max-n-content', type=float, help="Maxinum N content")
P_filter.add_argument('--max-rmask', type=float, help="Maxinum RepeatMask content")
P_filter.add_argument('--max-zero-depth', type=float, help="Maxinum Zero depth fraction")
P_filter.add_argument('--min-ref-coverage', type=float, help="Mininum reference coverage")
P_filter.add_argument('--max-ref-spread', type=float, help="Mininum reference spread")


P_sup = P_detect.add_argument_group("Additional Result")
P_sup.add_argument('--ideogram', help='USCS cytoBandIdeo.txt')
P_sup.add_argument('--plotcnv', action='store_true', help='Plot every chromosome and their CNV')
P_sup.add_argument('--annotation', action='store_true', help='Use ClassifyCNV to annotate CNV')



P_detect.set_defaults(func=_cmd_detect)

# version ---------------------------------------------------------------------

def print_version(_args):
    """Display this program's version."""
    print(__version__)


P_version = AP_subparsers.add_parser('version', help=print_version.__doc__)
P_version.set_defaults(func=print_version)


def mod_params(args):
    if args.min_gc:
        params.MIN_GC = args.min_gc
    if args.max_gc:
        params.MAX_GC = args.max_gc
    if args.min_alignability:
        params.MIN_ALIGNABILITY = args.min_alignability
    if args.max_n_content:
        params.MAX_N_CONTENT = args.max_n_content
    if args.max_rmask:
        params.MAX_RMASK = args.max_rmask
    if args.max_zero_depth:
        params.MAX_ZERO_DEPTH = args.max_zero_depth
    if args.min_ref_coverage:
        params.MIN_REF_COVERAGE = args.min_ref_coverage
    if args.max_ref_spread:
        params.MAX_REF_SPREAD = args.max_ref_spread



# _____________________________________________________________________________
# Shim for command-line execution

def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)

def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    args.func(args)

