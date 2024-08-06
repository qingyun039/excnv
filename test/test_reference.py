from cnvsq import reference
from cnvsq import coverage

genome = '/home/chenshulin/Project/Data/gatk-bundle-hg19/ucsc.hg19.fasta'
bigwig = ''

# ref = reference.do_reference_flat(10000)
# reference.probes2bed(ref)
# print(ref.data)

bamfiles = [ f'/data/oncology/yousheng/analysis/F350002939_20220902/cnvseq/align/YK05CS{i}.bam' for i in range(1, 6) ]
ref = reference.do_reference(bamfiles, 10000, fa_fname=genome)
print(ref.data)

# ref = reference.do_reference_flat(10000, fa_fname=genome)
# print(ref.data)

# ref = reference.do_reference_flat(10000, fa_fname=genome, bw_fname=bifwig)
# print(ref.data)
