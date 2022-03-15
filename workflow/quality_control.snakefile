## 1st step: Generate quality control statistics (based on fastqc and MultiQC)

from snakemake.utils import min_version

# Setup parameters
min_version("6.6.0")
# loading parameters from the configure file
SAMPLE_THREADS = config["sample_threads"]
FASTQDIR = config.get("fastqdir", "fastq")
FASTQPOSTFIX = config.get("fastqpostfix", "")

# fetching all samples
SAMPLE_LIST = config["metadata_file"]
with open(SAMPLE_LIST, "r") as samp_f:
    SAMPLES = [x.strip().split('\t')[0] for x in samp_f.readlines()[1:]]
    SAMPLES = list(x for x in SAMPLES if x)

rule all:
    input:
        "1.Quality_control/Sequence_QC/MultiQC_report.html"

##### Modules #####
include: "rules/quality_control.smk"
