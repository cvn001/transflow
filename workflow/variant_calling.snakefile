## 2nd step: variants calling (based on MTB pan-genome & bwa + GATK3)

import os
import re
import math
import subprocess
from snakemake.utils import min_version

min_version("6.6.0")

# Setup parameters
##### Config file and sample list file #####
# load configs
SAMPLE_THREADS = config.get("sample_threads", 1)
FASTQDIR = config.get("fastqdir", "fastq")
FASTQPOSTFIX = config.get("fastqpostfix", "")
# get sample list from sample_list.txt
SAMPLE_LIST = config.get("mtbc_samples", "2.MTBC_identification/MTBC_samples.txt")
with open(SAMPLE_LIST, "r") as samp_f:
    SAMPLES = [x.strip().split('\t')[0] for x in samp_f.readlines()]
    SAMPLES = list(x for x in SAMPLES if x)
# print('Total number of MTBC samples: {}'.format(len(SAMPLES)))
# print('MTBC samples: ', SAMPLES)

RESOURCE = srcdir('resources/')
SCRIPTS = srcdir("scripts/")
# Load Reference genome, default is a pre-built pan-genome
REF_GENOME = RESOURCE + 'reference/pangenomeMTB_consensus.fasta'
GENOME_FILE = config.get("genome_file", REF_GENOME)

GENOME_FDIR = os.path.dirname(GENOME_FILE)
GENOME = os.path.splitext(os.path.basename(GENOME_FILE))[0]

AF_threshold = config.get("allele_frequency_threshold", 0.75)
MQ_threshold = config.get("mapping_quality_threshold", 10)
DP_threshold = config.get("depth_threshold", 5)

# The prefix of output files.
OUTFILENAME = config.get("output_prefix", "samples")

##### Helper functions #####

def index_N_bases_from_fai(fai_file):
    with open(fai_file,'r') as f:
        first_line = f.readline()
        gl = int(first_line.split('\t')[1])
        gl = math.floor((math.log2(gl) / 2) - 1.5)
        gl = min(14, gl)
        return int(gl)


rule all:
    input:
        "1.Quality_control/Alignment_QC/MultiQC_report.html",
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_indels.vcf", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_snps.vcf", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_svs.vcf", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_uncovered.bed", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_uncalled.bed", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_deletions.bed", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_snps_amb-lowqual.bed", smp=SAMPLES, genome=GENOME),
        expand("3.SNP_calling/{smp}/{genome}/{smp}_vars_snps_lowqual.bed", smp=SAMPLES, genome=GENOME)  


##### Modules #####
include: "rules/mapping.smk"
include: "rules/bam_file_modification.smk"
include: "rules/alignment_qc.smk"
include: "rules/indices.smk"
include: "rules/snp_calling.smk"
include: "rules/vcf_filtering.smk"
