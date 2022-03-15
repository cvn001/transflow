library(yaml)

############ parse basic parameters ############ 
args <- commandArgs(trailingOnly = TRUE)
rmd_script <- args[1]
output <- args[2]
config_file <- args[3]
############ parse other parameters ############ 
config <- yaml.load_file(config_file)
config <- config[["dictitems"]]
kraken <- config$kraken
kraken_cutoff <- config$kraken_cutoff
multiqc <- config$multiqc
mapping_quality <- config$mapping_quality
allele_frequency <- config$allele_frequency
mapping_depth <- config$mapping_depth
snp_diff_matrix <- config$snp_diff_matrix
snp_diff_heatmap <- config$snp_diff_heatmap
snp_diff_histogram <- config$snp_diff_histogram
cluster_table <- config$cluster_table
clustering_visualization <- config$clustering_visualization
cluster_member_visualization <- config$cluster_member_visualization
show_network <- config$show_network
network_image <- config$network_image
risk_file <- config$risk_file
show_risk <- config$show_risk
cluster_method <- config$cluster_method
snp_cutoff <- config$snp_cutoff
trans_cutoff <- config$trans_cutoff

# run Rmarkdown script
rmarkdown::render(
  rmd_script,
  output_file = basename(output),
  output_dir = dirname(output),
  params = list(
    kraken = kraken,
    kraken_cutoff = kraken_cutoff,
    multiqc = multiqc,
    mapping_quality = mapping_quality,
    mapping_depth = mapping_depth,
    allele_frequency = allele_frequency,
    snp_diff_matrix = snp_diff_matrix,
    snp_diff_heatmap = snp_diff_heatmap,
    snp_diff_histogram = snp_diff_histogram,
    clustering_visualization = clustering_visualization,
    cluster_member_visualization = cluster_member_visualization,
    show_network = show_network,
    network_image = network_image,
    risk_file = risk_file,
    show_risk = show_risk,
    cluster_method = cluster_method,
    cluster_table = cluster_table,
    snp_cutoff = snp_cutoff,
    trans_cutoff = trans_cutoff
  )
)
