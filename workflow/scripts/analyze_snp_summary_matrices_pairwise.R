suppressWarnings(suppressMessages(require(argparse)))
suppressWarnings(suppressMessages(require(ggpubr)))
suppressWarnings(suppressMessages(library(GenomicRanges)))

# calculate hamming distance between columns in X
hamming <- function(X) {
  uniqs <- unique(as.vector(X))
  U <- X == uniqs[1]
  H <- t(U) %*% U
  for (uniq in uniqs[-1]) {
    U <- X == uniq
    H <- H + t(U) %*% U
  }
  nrow(X) - H
}

# calculate hamming distance between binary encoded columns in X
hamming_binary <- function(X) {
  D <- t(1 - X) %*% X
  D + t(D)
}

# get count of TRUE-TRUE pairs in binary encoded columns in X
true_binary <- function(X) {
  D <- t(X) %*% X
  diag(D) <- 0
  D
}

parser <- ArgumentParser(description='Calculate pairwise distance from data parsed from VCF and BED files.')
parser$add_argument('--input_data', dest='input_data', required='TRUE', help='R object generated with the get_snp_summary.R script.')
parser$add_argument('--output_dir', dest='output_dir', required='TRUE', help='Output directory for distance tables.')
parser$add_argument('--output_prefix', dest='output_prefix', required='TRUE', help='Prefix for distance tables.')
parser$add_argument('--exclude_regions', dest='exclude_bed', help='BED file with regions that should not be considered in distance calculation.')

# read input args
args <- parser$parse_args()
input_data <- args$input_data
output_dir <- args$output_dir
output_prefix <- args$output_prefix
additional_exclude_bed <- args$exclude_bed

# input data must be generated with "get_snp_summary.R"
cat("[1] Data loading.\n")
suppressWarnings(suppressMessages(load(input_data)))
genomelen <- as.numeric(genomelen)

# regions that should not be considered for distance (e.g. resistance mutations)
additional_exclude_pos <- c()
if (!is.null(additional_exclude_bed)) {
  bed <- tryCatch(read.table(additional_exclude_bed, stringsAsFactors = F, comment.char = "#", header=F), error = function(e) {data.frame()})
  if (ncol(bed) == 3) {
    colnames(bed) <- c("chr", "start", "end")
    bed[,1] <- genome
    bedasgranges <- makeGRangesFromDataFrame(bed, starts.in.df.are.0based = TRUE, seqinfo = si)
    additional_exclude_pos <- sapply(unique(unlist(mapply(start(bedasgranges), end(bedasgranges), FUN = ":"))), as.character)  
    exclude_pos <- intersect(additional_exclude_pos, colnames(snp_matrix))
    snp_matrix <- snp_matrix[,!(colnames(snp_matrix) %in% exclude_pos)]
    cat("[*] Excluded positions from matrix.\n")
  }
}

matrix_output <- paste0(output_dir, "/", output_prefix, "_snp_matrix.txt")
write.table(snp_matrix, matrix_output, row.names = T, col.names = NA, quote = F, sep = "\t")
cat("[2] Starting analysis. This may take a while...\n")


# hamming distance between columns including "X" for uncovered and substracting differences resulting from comparing A, C, G, T to "X"
pairwise_rel_matrix_m <- hamming(as.matrix(t(snp_matrix))) - hamming_binary(t(as.matrix(snp_matrix)) == "X")


if (length(additional_exclude_pos) > 0) {
  exclude_list <- lapply(exclude_list, setdiff, y = additional_exclude_pos)
}

# get all positions that are not covered in ALL samples
uncov_in_all <- setdiff(1:genomelen, additional_exclude_pos)
for (i in 1:length(exclude_list)) {
  uncov_in_all <- intersect(uncov_in_all, exclude_list[[i]])
}


# get all positions that are not covered in ANY samples
uncov_in_any <- c()
for (i in 1:length(exclude_list)) {
  uncov_in_any <- union(uncov_in_any, exclude_list[[i]])
}

# get all positions that are not covered in ANY samples, but not uncovered in all samples
# (uncovered positions that are covered in at least one sample)
uncov_cov_in_any <- setdiff(uncov_in_any, uncov_in_all)
uncov_cov_in_any <- setdiff(uncov_cov_in_any, additional_exclude_pos)
coveragem_length <- length(uncov_cov_in_any)

# create matrix for these positions, marking uncovered positions with "TRUE"
coverage_m <- matrix(F, nrow = coveragem_length, ncol = length(all_ids))
colnames(coverage_m) <- all_ids
rownames(coverage_m) <- sapply(uncov_cov_in_any, as.character)

for (id in all_ids) {
  pos <- exclude_list[[id]]
  pos <- intersect(pos, uncov_cov_in_any)
  pos <- sapply(pos, as.character)
  coverage_m[pos, id] <- T
}

# get count of TRUE - FALSE & FALSE - TRUE pairs with hamming distance,  and  FALSE - FALSE pairs with matrix multiplication
# to find all positions that are not covered when comparing two samples
# add number of positions that are uncovered in all samples
excluded_m <- hamming_binary(coverage_m) + true_binary(coverage_m) + length(uncov_in_all)
diag(excluded_m) <- 0

# divide differences by comparable regions (covered in both samples) and project unto genome length
coverage_m <- NULL
if (is.null(genomegaps)) {
  pairwise_rel_matrix_m_rel <- (pairwise_rel_matrix_m / (genomelen - length(additional_exclude_pos) - excluded_m)) * genomelen
} else {
  cat("[3] Integrating pan-genome gaps.\n")
  # if uncovered positions overlap between two samples and overlap with gaps in whole genome alignment of pangenome:
  # project differences unto shortened genome length 

  coverage_m_genomegaps <- intersect(uncov_cov_in_any, genomegaps) # get gaps uncovered in ANY sample
  uncov_all_genomgegaps_nr <- length(intersect(uncov_in_all, genomegaps)) # get gaps uncovered in ALL samples
  # create same matrix as before but only mark positions as uncovered if they overlap gaps in WGA
  coverage_m_genome <- matrix(F, nrow = length(coverage_m_genomegaps), ncol=length(all_ids))
  colnames(coverage_m_genome) <- all_ids
  rownames(coverage_m_genome) <- sapply(coverage_m_genomegaps, as.character)
  for (id in all_ids) { 
    pos <- exclude_list[[id]]
    pos <- intersect(pos, coverage_m_genomegaps)
    pos <- sapply(pos, as.character)
    coverage_m_genome[pos, id] <- T
  }
  # calculate shortened genome length for each comparison and use it for projecting difference count
  pairwise_intersect_genomegap_overlap <- genomelen - length(additional_exclude_pos) - (true_binary(coverage_m_genome) + uncov_all_genomgegaps_nr)
  diag(pairwise_intersect_genomegap_overlap) <- 0
  pairwise_rel_matrix_m_rel <- (pairwise_rel_matrix_m / (genomelen - length(additional_exclude_pos) - excluded_m)) * pairwise_intersect_genomegap_overlap
}
cat('[4] Building pairwise distance outputs.\n')
dir.create(file.path(output_dir), showWarnings = FALSE)

count_file <- paste0(output_dir, "/", output_prefix, "_pairwise_count.txt")
write.table(pairwise_rel_matrix_m, count_file, row.names = T, col.names = NA, quote = F, sep = "\t")
distance_file <- paste0(output_dir, "/", output_prefix, "_pairwise_distance_matrix.txt")
write.table(pairwise_rel_matrix_m_rel, distance_file, row.names = T, col.names = NA, quote = F, sep = "\t")
