#!/usr/bin/env Rscript
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(argparse)))

parser <- ArgumentParser(description = 'Parsing VCF files and BED files in preparation for pairwise distance calculation.')
parser$add_argument('--input_dir', dest = 'input_dir', required = 'TRUE', help = 'Path to SNP_Calling directory.')
parser$add_argument('--id_file', dest = 'id_f', required = 'TRUE', help = 'Sample ids to be analysed, one per line.')
parser$add_argument('--output_file', dest = 'output_f', required = 'TRUE', help = 'Full path and name of output Rdata file.')
parser$add_argument('--fastaidx_file', dest = 'faidx_f', required = 'TRUE', help = 'Fasta index file (.fasta.fai) of Fasta file used for alignment')
parser$add_argument('--genome_gaps_file', dest = 'genomegaps_f', required = 'TRUE', help = 'BED file with list of regions in pangenome with gaps in any genome of pangenome. Provide an empty file if a conventional reference genome was used.')

args <- parser$parse_args()

input_dir <- args$input_dir
id_f <- args$id_f
output_f <- args$output_f
faidx_f <- args$faidx_f
genomegaps_f <- args$genomegaps_f

genomeinfo <- read.table(faidx_f, header = F, sep = "\t", stringsAsFactors = F)
if (nrow(genomeinfo) > 1) {
   stop("Found more than one line in Fasta index file. Unfortunately PANPASCO does not work with more than one sequence per fasta file.")
}

genome <- sub(pattern = "\\.fai$", replacement = "", basename(faidx_f))
genome <- sub(pattern = "\\.fasta$", replacement = "", genome)
genomelen <- as.numeric(genomeinfo[1,2])


info <- file.info(genomegaps_f)
if (info[1, "size"] == 0) {
  genomegaps_f <- NULL
}


# create seqinfo object for genomic ranges 
si <- Seqinfo(seqnames = c(genome), seqlengths = c(genomelen))
snp_matrix <- data.frame()

ref <- data.frame()
all_ids <- unlist(read.table(id_f, header = F, stringsAsFactors = F))

# read all vcf positions with bases and store bases in snp matrix
# also keep track of reference bases for each position
cat("[1] Reading VCF files.\n")
for (id in all_ids) {
  cat(".")
  vcf_file <- paste0(input_dir, "/", id, "/", genome, "/", id, "_vars_snps.vcf")
  vcf <- read.table(vcf_file, stringsAsFactors = F, comment.char = "#", header = F)
  for (i in 1:nrow(vcf)) {
    pos <- as.character(vcf[i,2])
    ref["ref", pos] <- vcf[i,4]
    snp_matrix[id, pos] <- vcf[i,5]
  }
}
cat('\n')
snp_positions <- colnames(snp_matrix)

# read all low quality regions from BED - files
# use Genomic Ranges for fast handling of intervals
cat("[2] Reading BED files.\n")
bedpostfixes <- c("uncalled", "uncovered", "vars_deletions", "vars_snps_amb-lowqual", "vars_snps_lowqual")
exclude_list <- list()
index <- 1
for (id in all_ids) {
  cat(index)
  exclude_list_cur <- GRanges()
  for (bedpostfix in bedpostfixes) {
    cat(".")
    bed_file <- paste0(input_dir, "/", id, "/", genome, "/", id, "_", bedpostfix, ".bed")
    bed_data <- read.table(bed_file, stringsAsFactors = F, comment.char = "#", header = F)
    bed <- tryCatch(bed_data, error = function(e) {data.frame()})
    if (ncol(bed) == 3) {
      colnames(bed) <- c("chr", "start", "end")
      bed[,1] <- genome
      bedasgranges <- makeGRangesFromDataFrame(bed, starts.in.df.are.0based = TRUE, seqinfo = si)
      exclude_list_cur <- union(exclude_list_cur, bedasgranges)
    }
  }
  exclude_list[[id]] <- unique(unlist(mapply(start(exclude_list_cur), end(exclude_list_cur), FUN = ":")))
  index <- index + 1
}
cat('\n')
cat("[3] Adding low quality information to matrix.\n")
# insert 'X' into snp matrix for all low quality regions 
for (id in all_ids) {
  cat(".")
  excl_positions <- exclude_list[[id]]
  positions <- intersect(excl_positions, snp_positions)
  positions <- sapply(positions, as.character)
  snp_matrix[id, positions] <- "X"
}
cat('\n')
# insert reference base for positions without alt (from vcf file) or evidence of low quality (from bed files)
cat("[4] Inserting reference base for positions without alt (from vcf file) or evidence of low quality (from bed files).\n")
for (id in all_ids) {
  cat(".")
  na_pos <- colnames(snp_matrix[id, is.na(snp_matrix[id,])])
  snp_matrix[id, na_pos] <- ref["ref", na_pos]
}
cat('\n')
cat("[5] Reading intervals with gaps in whole genome alignment for the pangenome.\n")
# read intervals with gaps in whole genome alignment for pangenome
genomegaps <- NULL
if (!is.null(genomegaps_f)) {
  cat("[6] Integrating pan-genome gaps.\n")
  gap_data <- read.table(genomegaps_f, stringsAsFactors = F, comment.char = "#", header = F)
  genomegapsbed <- tryCatch(gap_data, error = function(e){NULL})
  if (!is.null(genomegapsbed)) {
    colnames(genomegapsbed) <- c("chr", "start", "end")
    genomegapsbed[,1] <- genome
    genomegaps <- disjoin(makeGRangesFromDataFrame(genomegapsbed, starts.in.df.are.0based = TRUE, seqinfo = si))
    genomegaps <- unlist(mapply(start(genomegaps), end(genomegaps), FUN = ":"))
  }
}
cat("[7] Saving SNP summary.\n")
save("snp_matrix", "exclude_list", "all_ids", "genome", "genomelen", "genomegaps", "ref", "si", file = output_f)

