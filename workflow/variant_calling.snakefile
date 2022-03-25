## 3rd step: variants calling (based on MTB pan-genome & bwa + GATK3)

import os
import re
import math
import subprocess
from snakemake.utils import min_version

# Setup parameters
min_version("6.6.0")

##### Config file and sample list file #####
# load configs

SAMPLE_THREADS = config["sample_threads"]
FASTQDIR = config.get("fastqdir", "fastq")
FASTQPOSTFIX = config.get("fastqpostfix", "")
# get sample list from sample_list.txt
SAMPLE_LIST = config.get("mtbc_samples", "2.MTBC_identification/MTBC_samples.txt")
with open(SAMPLE_LIST, "r") as samp_f:
    SAMPLES = [x.strip().split('\t')[0] for x in samp_f.readlines()]
    SAMPLES = list(x for x in SAMPLES if x)

FLASH_OVERLAP = config.get("flash_overlap", 10)
READ_MIN_LEN = config.get("trimmomatic_read_minimum_length", 50)
# Load Reference genome, default is a pre-built pan-genome
RESOURCE = srcdir('resources/')
REF_GENOME = RESOURCE + 'reference/pangenomeMTB_consensus.fasta'
GENOME_FILE = config.get("genome_file", REF_GENOME)
GENOME_FDIR = os.path.dirname(GENOME_FILE)
GENOME = os.path.splitext(os.path.basename(GENOME_FILE))[0]
ADATER = RESOURCE + 'trimmomatic_adapter/NexteraPE-PE.fa'
ADAPTER_FILE = config.get("trimmomatic_adapter_file", ADATER)
MIN_QUAL = config.get("trimmomatic_qual_slidingwindow", 15)

AF_threshold = config.get("allele_frequency_threshold", 0.75)
MQ_threshold = config.get("mapping_quality_threshold", 10)
DP_threshold = config.get("depth_threshold", 5)

CLIPPING = FLASH_OVERLAP / 2
SCRIPTS = srcdir("scripts/")
# The prefix of output files.
OUTFILENAME = config.get("output_prefix", "MTBC_samples")


##### Helper functions #####

def read_length_from_file(wildcards):
    # get average read length
    shell_cmd = " | awk 'BEGIN {ml=0} {if(NR%4==2) { x=length($1); if (x+0>ml+0) ml=x}} END {print ml}'"
    cmd_1 = "gzip -dc temp/adapt_clip/" + wildcards['smp'] + "/" + wildcards['smp'] + "_1.fastq.gz" + shell_cmd
    cmd_2 = "gzip -dc temp/adapt_clip/" + wildcards['smp'] + "/" + wildcards['smp'] + "_2.fastq.gz" + shell_cmd
    rl1 = int(subprocess.check_output(cmd_1, shell=True))
    rl2 = int(subprocess.check_output(cmd_2, shell=True))
    return max(rl1, rl2)


def read_length_from_name(wildcards):
    m = re.search('(\d+)bp',wildcards['smp'])
    if m is not None:
        rl = int(m.group(1))
    else:
        rl = read_length_from_file(wildcards)
    rl_clip = rl - CLIPPING
    return int(rl), int(rl_clip)


def read_length_from_name_70(wildcards):
    rl = read_length_from_name(wildcards)
    rl_70 = math.floor(rl[0] * 0.7)
    return int(rl_70)


def index_N_bases_from_fai(fai_file):
    with open(fai_file,'r') as f:
        first_line = f.readline()
        gl = int(first_line.split('\t')[1])
        gl = math.floor((math.log2(gl) / 2) - 1.5)
        gl = min(14, gl)
        return int(gl)


rule all:
    input:
        expand("temp/ready_forvar/{smp}/{genome}/{smp}_ready.bam", smp=SAMPLES, genome=GENOME),
        expand("1.Quality_control/Alignment_QC/MultiQC_report.html"),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_indels.vcf", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_snps.vcf", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_svs.vcf", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_uncovered.bed", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_uncalled.bed", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_deletions.bed", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_snps_amb-lowqual.bed", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_snps_lowqual.bed", smp=SAMPLES, genome=GENOME)  


##### Modules #####
include: "rules/clipping.smk"
include: "rules/transition.smk"
include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/bam_file_modification.smk"
include: "rules/alignment_qc.smk"
include: "rules/indices.smk"
include: "rules/snp_calling.smk"
include: "rules/vcf_filtering.smk"
