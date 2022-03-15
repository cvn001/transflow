suppressWarnings(suppressMessages(library(pheatmap)))
suppressWarnings(suppressMessages(library(ggpubr)))

cat("[1] Drawing pairwise snp distance matrix heatmap.\n")
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]
show_colnames_arg <- FALSE
if (args[3] %in% c("TRUE", "True", "true", "T", "t")) {
  show_colnames_arg <- TRUE
}
show_rownames_arg <- FALSE
if (args[4] %in% c("TRUE", "True", "true", "T", "t")) {
  show_rownames_arg <- TRUE
}
prefix <- args[5]
heatmap_name <- paste0(prefix, "_pairwise_distance_heatmap")
output_heatmap_pdf <- file.path(output_dir, paste0(heatmap_name, ".pdf"))
output_heatmap_png <- file.path(output_dir, paste0(heatmap_name, ".png"))
snp_diff <- read.table(input_file, header = T, check.names = F, sep = "\t", row.names = 1) # nolint
cat("[*] Exporting heatmap image.\n")
# export to a pdf file
pheatmap(snp_diff, show_rownames = show_rownames_arg, fontsize = 3.5, border_color = NA,
  show_colnames = show_colnames_arg, filename = output_heatmap_pdf)
# export to a png file
pheatmap(snp_diff, show_rownames = show_rownames_arg, fontsize = 3.5, border_color = NA,
         show_colnames = show_colnames_arg, filename = output_heatmap_png)
cat("[2] Drawing pairwise snp distance histogram plot.\n")
snp_diff <- as.matrix(snp_diff)
snp_diff[lower.tri(snp_diff, diag = T)] <- NA
long_df <- na.omit(data.frame(as.table(snp_diff)))
names(long_df) <- c("A", "B", "Distance")
line_pos <- 50
long_df$Group <- ifelse(long_df$Distance > line_pos, 'far', 'close')
# draw histogram
p <- gghistogram(long_df, x = "Distance", rug = T, bins = 90, ylab = "Count", fill = "white", color = "Group", palette = c("red", "#0F425CFF")) + # nolint
  geom_vline(xintercept = line_pos, col = "red", linetype = "dashed") + xlab("Pairwise SNP distance") + ylab("Frequency") + # nolint
  theme(legend.position = "none")
histogram_name <- paste0(prefix, "_pairwise_distance_distribution")
cat("[*] Exporting histogram image.\n")
output_histogram_pdf <- file.path(output_dir, paste0(histogram_name, ".pdf"))
output_histogram_png <- file.path(output_dir, paste0(histogram_name, ".png"))
ggsave(output_histogram_pdf, width = 7, height = 5)
ggsave(output_histogram_png, width = 7, height = 5, dpi = 300)
