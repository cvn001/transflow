## 4th step: generating a summary report (html format)

import os
from snakemake.utils import min_version

min_version("6.6.0")

# Setup parameters
SCRIPTS = srcdir("scripts/")
WORK_DIR = os.getcwd()
# Report parameters
KRAKEN_CUTOFF = config.get("kraken_cutoff", 90)
AF_threshold = config.get("allele_frequency_threshold", 0.75)
MQ_threshold = config.get("mapping_quality_threshold", 20)
DP_threshold = config.get("depth_threshold", 10)
SNP_CUTOFF = config.get("snp_threshold", 12)
TRANS_CUTOFF = config.get("transmission_threshold", 15)
LAMBDA = config.get("clock_rate", 1.0)
BETA = config.get("transmission_rate", 2.0)
MTBC_READS_ONLY = config.get("MTBC_reads_only", True)
OUTFILENAME = config.get("output_prefix", "samples")
# Transmission clustering method
METHOD = config.get("method", "trans")
if METHOD in ['snp', 'SNP']:
    METHOD = 'SNP'
    METHOD_DIR = 'SNP_based_method'
    CUTOFF = SNP_CUTOFF
else:
    METHOD = 'trans'
    METHOD_DIR = 'Transmission_based_method'
    CUTOFF = TRANS_CUTOFF

rule all:
    input:
        expand("7.Summary_report/{output}/summary_report.html", output=METHOD_DIR)

##### Modules #####
include: "rules/summary_report.smk"
