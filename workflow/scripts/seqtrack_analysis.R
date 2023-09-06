suppressWarnings(suppressMessages(require(argparse)))
suppressWarnings(suppressMessages(require(GGally)))
suppressWarnings(suppressMessages(require(ggplot2)))
suppressWarnings(suppressMessages(require(adegenet)))
suppressWarnings(suppressMessages(require(gridExtra)))
suppressWarnings(suppressMessages(require(lubridate)))

rm(list = ls())
parser <- ArgumentParser(description = 'Using SeqTrack to infer the transmission network.')
parser$add_argument('-i', '--distance', dest = 'distance', required = 'TRUE', help = 'Input symmetric distance matrix to use for the detection.')
parser$add_argument('-s', '--sample', dest = 'sample', required = 'TRUE', help = 'Input sample list to use for the detection.')
parser$add_argument('-m', '--meta', dest = 'meta', required = 'TRUE', help = 'Metadata of samples.')
parser$add_argument('-c', '--coord', dest = 'coord', required = 'FALSE', default = 'false', help = 'Using latitude and longitude.')
parser$add_argument('-o', '--output_dir', dest = 'output_dir', required = 'TRUE', help = 'Output directory.')
# 
## read input arguments
args <- parser$parse_args()
input_distance <- args$distance
input_sample <- args$sample
output_dir <- args$output_dir
input_meta <- args$meta
use_coord <- args$coord

set.seed(12345)

meta_data <- read.table(input_meta, header = T, sep = '\t', stringsAsFactors = F, check.names = F)
meta_data$sample <- as.character(meta_data$sample)
samples <- read.csv(input_sample, header = F, check.names = F)
colnames(samples) <- c('sample')
samples$sample <- as.character(samples$sample)
distance <- read.table(input_distance, sep = '\t', header = T, row.names = 1, stringsAsFactors = F, check.names = F)
distance <- round(as.matrix(distance), 1)
if (use_coord %in% c('true', 'True', 'TRUE', 't', 'T')) {
  use_coord <- T
  cols <- c('sample', 'date', 'latitude', 'longitude')
  col_string <- paste(cols, collapse = ', ')
  check <- ifelse(cols %in% colnames(meta_data), 1, 0)
  if (0 %in% check) {
    print(paste0("Error: not valid columns [", col_string, "] found in metadata file."))
    quit(status = 1)
  } else {
    use_prox <- T
  }
} else {
   use_coord <- F
  cols <- c('sample', 'date')
  col_string <- paste(cols, collapse = ', ')
  check <- ifelse(cols %in% colnames(meta_data), 1, 0)
  if (0 %in% check) {
    print(paste0("Error: not valid columns [", col_string, "] found in metadata file."))
    quit(status = 1)
  } else {
    use_prox <- F
  }
}
meta_data <- subset(meta_data, meta_data$sample %in% samples$sample)
new_order <- sapply(samples$sample, function(x, df) {which(df$sample == x)}, df = meta_data)
meta_data <- meta_data[new_order,]
D <- distance[samples$sample, samples$sample]
# Convert four digit year values to class Date, "2007" to "2007-01-01"

# dates <- as.Date(ISOdate(meta_data$date, 1, 1))
dates <- parse_date_time(meta_data$date, c('ymd', 'ym', 'y'))
distmat <- as.matrix(D)
cases <- meta_data$sample

# run seqtrack function
if (use_prox) {
  xy <- cbind(meta_data$longitude, meta_data$latitude)
  temp <- as.matrix(dist(xy))
  M <- 1 * (temp < 1e-10)
  res <- seqTrack(distmat, x.names = cases, x.dates = dates, prox.mat = M)
} else {
  res <- seqTrack(distmat, x.names = cases, x.dates = dates)
}
new_res <- data.frame('sample' = rownames(res), res)
annotate_table <- new_res[, c('id', 'sample', 'date')]

width <- c(2, 1)
size <- 5
if (use_prox) {
  meta_data <- meta_data[, c('sample', 'latitude', 'longitude')]
  annotate_table <- merge(annotate_table, meta_data, by = 'sample')
  annotate_table <- annotate_table[, c('id', 'sample', 'date', 'latitude', 'longitude')]
  width <- c(2, 1)
  size <- 4
}
annotate_table <- annotate_table[order(annotate_table$id), , drop = F]
res.file <- paste0(output_dir, '/', 'seqtrack_result.txt')
write.table(data.frame('sample' = rownames(res), res), file = res.file, sep = '\t', quote = F, row.names = F, na = '')
from.old <- res$ances
to.old <- res$id
weight <- res$weight
isNotNA <- !is.na(from.old) & !is.na(to.old)
vnames <- sort(unique(c(from.old, to.old)))
from <- match(from.old, vnames)
to <- match(to.old, vnames)
dat <- data.frame(from, to, weight, stringsAsFactors = F)[isNotNA, , drop = F]
cyto_dir <- file.path(output_dir, 'for_cytoscape')
# export files for cytoscape visualization
dir.create(cyto_dir, showWarnings = F)
network_file <- file.path(cyto_dir, 'network.txt')
write.table(dat, file = network_file, sep = '\t', quote = F, row.names = F)
res$sample <- row.names(res)

if (use_coord) {
  node <- merge(res[, c('sample', 'id', 'date')], meta_data[, c('sample', 'longitude', 'latitude')])
} else {
  node <- res[, c('sample', 'id', 'date')]
}
node_file <- file.path(cyto_dir, 'node.txt')
write.table(node, file = node_file, sep = '\t', quote = F, row.names = F)

# Draw transmission network using ggnet2
dat$weight <- round(dat$weight, digits = 0)
net <- try(network::network(dat), silent = TRUE)
if (network::is.network(net)) {
  p <- ggnet2(dat, label = T, arrow.size = 5, arrow.gap = 0.05, size = 0,
        label.size = 4, color = 'steelblue', mode = 'kamadakawai',
        edge.color = 'grey30', edge.label = 'weight', edge.label.size = 2.5) +
    guides(size = "none") +
    geom_point(aes(color = color), size = 10, color = "white") +
    geom_point(aes(color = color), size = 10, alpha = 0.5) +
    geom_point(aes(color = color), size = 8) +
    geom_text(aes(label = label), color = "white", fontface = "bold")
  t <- tableGrob(annotate_table, rows = NULL, theme = ttheme_default(base_size = size))
  g <- arrangeGrob(p, t, nrow = 1, widths = width)
  pdf_file <- file.path(output_dir, 'transmission_network.pdf')
  png_file <- file.path(output_dir, 'transmission_network.png')
  suppressMessages(ggsave(file = pdf_file, g, width = 6, height = 4))
  suppressMessages(ggsave(png_file, g, dpi = 500, width = 6, height = 4))
} else {
  cat("[*] WARNING, could not coerce net to a network object, skipping.\n")
}
