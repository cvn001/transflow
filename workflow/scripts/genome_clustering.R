suppressWarnings(suppressMessages(require(argparse)))
suppressWarnings(suppressMessages(library(transcluster)))
suppressWarnings(suppressMessages(library(stats)))
suppressWarnings(suppressMessages(library(lubridate)))

options(warn = -1)
rm(list = ls())

parser <- ArgumentParser(description = 'Using transcluster package to identify outbreak clusters using inferred transmissions.')
parser$add_argument('--distance', dest = 'distance', required = 'TRUE', help = 'Pairwise distance matrix.')
parser$add_argument('--date', dest = 'date', required = 'FALSE', help = 'Sample collection date.')
parser$add_argument('--lambda', dest = 'lambda', default = 1.0, help = '(optional) Clock rate [default 1.0].')
parser$add_argument('--beta', dest = 'beta', default = 2.0, help = '(optional) Transmission rate [default 2.0].')
parser$add_argument('--output_dir', dest = 'output_dir', required = 'TRUE', help = 'Output directory for clustering results.')
parser$add_argument('--prefix', dest = 'prefix', required = 'TRUE', help = 'Prefix for distance tables.')
parser$add_argument('--decimal_date', dest = 'decimal_date', required = 'TRUE', help = 'Convert date to decimal format.')
parser$add_argument('--method', dest = 'method', default = 'transmission', help = '(optional) Transmission clustering method [default transmission].')
parser$add_argument('--cutoff', dest = 'cutoff', default = 12, help = '(optional) SNP or Transmission threshold [default 12].')

# read input args
args <- parser$parse_args()
distance_file <- args$distance
date_file <- args$date
my_lambda <- as.numeric(args$lambda)
my_beta <- as.numeric(args$beta)
output_dir <- args$output_dir
prefix <- args$prefix
decimal_date <- args$decimal_date
detection_method <- args$method
cutoff <- as.numeric(args$cutoff)

set.seed(12345)
# Loading pairwise SNP distance matrix data.
snp_df <- read.csv(distance_file, sep = '\t', header = T, row.names = 1, check.names = F)

if (!file.exists(date_file)) {
  print(paste0("Error: sample collection date file does not exist.\n"))
  quit(status = 1)
}
df <- read.csv(date_file, sep = '\t', header = T, check.names = F, comment.char = "#")
cols <- c('sample', 'date')
col_string <- paste(cols, collapse = ', ')
check <- ifelse(cols %in% colnames(df), 1, 0)
if (0 %in% check) {
  print(paste0("Error: not valid columns [", col_string, "] found in metadata file.\n"))
  quit(status = 1)
}
df <- df[, c('sample', 'date')]
rownames(df) <- df$sample
if (decimal_date %in% c('TRUE', 'True', 'true', 'T', 't')) {
  df$ds <- ymd(df$date)
  df$date <- decimal_date(df$ds)
}
inter_ids <- intersect(row.names(snp_df), df$sample)
df <- df[inter_ids,]
snp_df <- snp_df[inter_ids, inter_ids]
ids <- df$sample
dates <- df$date
snpMatrix <- round(as.matrix(snp_df))
cat("[1] Loading pairwise SNP distance matrix and date data.\n")
myModel <- createModel(ids, dates, snpMatrix)
dir.create(output_dir, showWarnings = F)
setwd(output_dir)
if (detection_method == 'SNP') {
  cat(paste0("[2] Setting the SNP threshold: [", cutoff, "].\n"))
  myModel <- setSNPThresholds(myModel, thSNP = cutoff)
  mySNPClusters <- makeSNPClusters(myModel, nameBase = prefix, writeFile = T)
} else {
  cat("[*] Beyond the SNP threshold: identifying outbreak clusters using inferred transmissions.\n")
  myModel <- setParams(myModel, lambda = my_lambda, beta = my_beta)
  cat(paste0("[2] Setting the transmission threshold: [", cutoff, "].\n"))
  myModel <- setTransThresholds(myModel, thTrans = cutoff)
  myTransClusters <- makeTransClusters(myModel, nameBase = prefix, writeFile = T)
}
