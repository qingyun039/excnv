"""Hard-coded parameters for CNVkit. These should not change between runs."""
# Filter thresholds used in constructing the reference (log2 scale)
MIN_REF_COVERAGE = -5.0
MAX_REF_SPREAD = 0.8  # 可以据此过滤多态性区域
NULL_LOG2_COVERAGE = -20.0

# GC含量对WGS影响不大，主要影响捕获步骤，原始是0.3 ~ 0.7
MIN_GC = 0.3
MAX_GC = 0.7

# 使用BWA比对时多重比对随机选择一个比对位置, 如要设置，可设为0.5
MIN_ALIGNABILITY = 0.6
MAX_N_CONTENT = 0.8
MAX_RMASK = 0.9
MAX_ZERO_DEPTH = 0.1
# Fragment size for paired-end reads
INSERT_SIZE = 250

# Target/bin names that are not meaningful gene names
# (In some UCSF panels, "CGH" probes denote selected intergenic regions)
IGNORE_GENE_NAMES = ("-", ".", "CGH")
ANTITARGET_NAME = "Antitarget"
ANTITARGET_ALIASES = (ANTITARGET_NAME, "Background")

CHROMS = [ 'chr'+str(i) for i in range(1,23) ] + ['chrX', 'chrY']
