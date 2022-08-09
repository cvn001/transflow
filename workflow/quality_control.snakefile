## 1st step: Generate quality control statistics (based on fastqc and MultiQC)

from snakemake.utils import min_version

min_version("6.6.0")

# Setup parameters
# loading parameters from the configure file
SAMPLE_THREADS = config.get("sample_threads", 1)
FASTQDIR = config.get("fastqdir", "fastq")
FASTQPOSTFIX = config.get("fastqpostfix", "")

RESOURCE = srcdir("resources/")
SCRIPTS = srcdir("scripts/")

REF_GENOME = RESOURCE + "reference/pangenomeMTB_consensus.fasta"
GENOME_FILE = config.get("genome_file", REF_GENOME)
GENOME_FDIR = os.path.dirname(GENOME_FILE)
GENOME = os.path.splitext(os.path.basename(GENOME_FILE))[0]

# ADATER = RESOURCE + 'trimmomatic_adapter/NexteraPE-PE.fa'
# ADAPTER_FILE = config.get("trimmomatic_adapter_file", ADATER)

# Kraken v1 parameters
KRAKEN_FILTER = config.get("kraken_cutoff", 50)
MTBC_READS_ONLY = config.get("MTBC_reads_only", True)
if MTBC_READS_ONLY:
    print('[*] Only MTBC reads will be kept for following steps\n')
    
# Trimmomatic parameters
MIN_QUAL = config.get("trimmomatic_qual_slidingwindow", 15)
READ_MIN_LEN = config.get("trimmomatic_read_minimum_length", 50)

# Flash parameters
FLASH_OVERLAP = config.get("flash_overlap", 10)
CLIPPING = FLASH_OVERLAP / 2

# fetching all input samples
GLOB_FILES = config.get("glob_files", False)
if GLOB_FILES:
    print("[*] Parsing all samples by glob patterns.\n")
    SAMPLES = glob_wildcards(FASTQDIR + "/{sample}_" + FASTQPOSTFIX + "1.fastq.gz").sample
else:
    print("[*] Parsing all samples from the metadata file.\n")
    METADATA_File = config["metadata_file"]
    with open(METADATA_File, "r") as samp_f:
        SAMPLES = [x.strip().split('\t')[0] for x in samp_f.readlines()[1:]]
        SAMPLES = list(x for x in SAMPLES if x)
print("=> Total number of samples: {}\n".format(len(SAMPLES)))
print("=> Samples are: \n\n", SAMPLES)

##### Helper functions #####

def read_length_from_file(wildcards, MTBC_READS_ONLY):
    # get average read length
    shell_cmd = " | awk 'BEGIN {ml=0} {if(NR%4==2) { x=length($1); if (x+0>ml+0) ml=x}} END {print ml}'"
    if MTBC_READS_ONLY:
        ap = '_mtbc'
    else:
        ap = ''
    cmd_1 = "gzip -dc temp/adapt_clip/" + wildcards['smp'] + "/" + wildcards['smp'] + "_cleaning{}_1.fastq.gz".format(ap) + shell_cmd
    cmd_2 = "gzip -dc temp/adapt_clip/" + wildcards['smp'] + "/" + wildcards['smp'] + "_cleaning{}_2.fastq.gz".format(ap) + shell_cmd
    rl1 = int(subprocess.check_output(cmd_1, shell=True))
    rl2 = int(subprocess.check_output(cmd_2, shell=True))
    return max(rl1, rl2)


def read_length_from_name(wildcards, MTBC_READS_ONLY):
    m = re.search('(\d+)bp',wildcards['smp'])
    if m is not None:
        rl = int(m.group(1))
    else:
        rl = read_length_from_file(wildcards, MTBC_READS_ONLY)
    rl_clip = rl - CLIPPING
    return int(rl), int(rl_clip)


def read_length_from_name_70(wildcards, MTBC_READS_ONLY):
    rl = read_length_from_name(wildcards, MTBC_READS_ONLY)
    rl_70 = math.floor(rl[0] * 0.7)
    return int(rl_70)


rule all:
    input:
        "1.Quality_control/Sequence_QC/After_trimming/MultiQC_report.html",
        "1.Quality_control/Sequence_QC/Before_trimming/MultiQC_report.html",
        "2.MTBC_identification/MTBC_samples.txt"

##### Modules #####
include: "rules/adapter_clipping.smk"
include: "rules/transition.smk"
include: "rules/trimming.smk"
include: "rules/quality_control.smk"
include: "rules/mtbc_filtering.smk"
