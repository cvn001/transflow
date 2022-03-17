## 4th step: sample clustering and transmission events predicting

import os
from snakemake.utils import min_version

# Setup parameters
min_version("6.6.0")
# loading parameters from the configure file
META_DATA_FILE = config["metadata_file"]
SAMPLE_LIST = config.get("mtbc_samples", "2.MTBC_identification/MTBC_samples.txt")
with open(SAMPLE_LIST, "r") as samp_f:
    SAMPLES = [x.strip().split('\t')[0] for x in samp_f.readlines()]
    SAMPLES = list(x for x in SAMPLES if x)
WORK_DIR = os.getcwd()
OUTFILENAME = "samples"
SCRIPTS = srcdir("scripts/")
# Reference genome data
RESOURCE = srcdir('resources/')
REF_GENOME = RESOURCE + 'reference/MTB_pangenome_consensus.fasta'
GENOME_FILE = config.get("genome_file", REF_GENOME)
GENOME = os.path.splitext(os.path.basename(GENOME_FILE))[0]
REF_GENOME_BED = RESOURCE + 'reference/MTB_pangenome_genome_gaps.bed'
EXCLUDE_BED = RESOURCE + 'exclusion/MTB_pangenome_excluded_loci.bed'
GENOMEGAPSFILE = config.get("genomegaps_file", REF_GENOME_BED)
EXCLUDEREGIONS = config.get("exclude_regions_file", EXCLUDE_BED)
# Alignment statistics
COVERAGE_X1 = config.get("1x_coverage", 88)
COVERAGE_X10 = config.get("10x_coverage", 75)
# SNP parameters
VALIDITY = config.get("check_transmission_validity", "False")
SNP_CUTOFF = config.get("snp_threshold", 12)
# Transmission parameters
TRANS_CUTOFF = config.get("transmission_threshold", 15)
LAMBDA = config.get("clock_rate", 1.0)
BETA = config.get("transmission_rate", 2.0)
DECIMAL_DATE = config.get("decimal_date", "True")
COORDINATE = config.get("coordinate", "False")
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
NETWORK = config.get("network_inference", "True")
# Epidemiological characteristics
CHARACTERS = config.get("characteristics", "")

def check_exclude(wildcards):
    if EXCLUDEREGIONS:
        if not os.path.exists(EXCLUDEREGIONS):
            raise Exception('Missing input files for rule distance:\n{}'.format(EXCLUDEREGIONS))
        print('Exclude regions file: {}'.format(EXCLUDEREGIONS))
        return " ".join([" --exclude_regions", EXCLUDEREGIONS])
    return ""

basic_output = list()
basic_output.append(expand("4.SNP_distance/{outfilename}.RSave", outfilename=OUTFILENAME))
basic_output.append(expand("4.SNP_distance/{outfilename}_pairwise_distance_matrix.txt", outfilename=OUTFILENAME))
basic_output.append(expand("4.SNP_distance/{outfilename}_snp_matrix.txt", outfilename=OUTFILENAME))
basic_output.append(expand("4.SNP_distance/{outfilename}_pairwise_distance_distribution.pdf", outfilename=OUTFILENAME))
basic_output.append(expand("4.SNP_distance/{outfilename}_pairwise_distance_distribution.png", outfilename=OUTFILENAME))
basic_output.append(expand("4.SNP_distance/{outfilename}_pairwise_distance_heatmap.pdf", outfilename=OUTFILENAME))
basic_output.append(expand("4.SNP_distance/{outfilename}_pairwise_distance_heatmap.png", outfilename=OUTFILENAME))
basic_output.append(expand("5.Transmission_cluster/{output}/transmission_clusters.txt", output=METHOD_DIR))
basic_output.append(expand("5.Transmission_cluster/{output}/transmission_detection_run_info.txt", output=METHOD_DIR))
basic_output.append(expand("5.Transmission_cluster/{output}/clustering_visualization.pdf", output=METHOD_DIR))
basic_output.append(expand("5.Transmission_cluster/{output}/clustering_visualization.png", output=METHOD_DIR))
basic_output.append(expand("5.Transmission_cluster/{output}/cluster_member_visualization.pdf", output=METHOD_DIR))
basic_output.append(expand("5.Transmission_cluster/{output}/cluster_member_visualization.png", output=METHOD_DIR))

if CHARACTERS:
    basic_output.append(expand("6.Risk_factor/{output}/univariable_regression_result.pdf", output=METHOD_DIR))
    basic_output.append(expand("6.Risk_factor/{output}/univariable_regression_result.png", output=METHOD_DIR))
    basic_output.append(expand("6.Risk_factor/{output}/characteristic_data.txt", output=METHOD_DIR))

rule all:
    input:
        basic_output

##### Modules #####
include: "rules/pairwise_distance.smk"
include: "rules/transmission_detection.smk"
include: "rules/risk_factor.smk"
