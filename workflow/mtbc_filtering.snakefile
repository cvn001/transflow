## 2nd step: MTBC filtering (based on Kraken)

from snakemake.utils import min_version

# Setup parameters
min_version("6.6.0")
# loading parameters from the configure file
SAMPLE_THREADS = config["sample_threads"]
FASTQDIR = config.get("fastqdir", "fastq")
FASTQPOSTFIX = config.get("fastqpostfix", "")
KRAKEN_FILTER = config.get("kraken_cutoff", 90)
SCRIPTS = srcdir("scripts/")
# fetching all samples
SAMPLE_LIST = config["metadata_file"]
with open(SAMPLE_LIST, "r") as samp_f:
    SAMPLES = [x.strip().split('\t')[0] for x in samp_f.readlines()[1:]]
    SAMPLES = list(x for x in SAMPLES if x)

rule all:
    input:
        "2.MTBC_identification/MTBC_samples.txt"

##### Modules #####
include: "rules/mtbc_filtering.smk"
