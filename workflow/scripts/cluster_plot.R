suppressWarnings(suppressMessages(require(argparse)))
suppressWarnings(suppressMessages(require(ggpubr)))
suppressWarnings(suppressMessages(require(ggsci)))

options(warn = -1)
rm(list = ls())

cat("[0] Drawing visualizations of transmission clustering results.\n")
parser <- ArgumentParser(description = 'Using transcluster package to identify outbreak clusters using inferred transmissions.')
parser$add_argument('--working_dir', dest = 'working_dir', required = 'TRUE', help = 'Working directory for clustering results.')
parser$add_argument('--input', dest = 'input', required = 'TRUE', help = 'Input file.')

####### read input args
args <- parser$parse_args()
working_dir <- args$working_dir
input_file <- args$input
####### 1st plot: a pie chart of the transmission clustering results
cat("[1] First pie chart...\n")
df <- read.table(input_file, sep = '\t', header = T, row.names = 1)
cluster_df <- as.data.frame(table(df[,2]))
colnames(cluster_df) <- c('type', 'count')
t_title <- 'Visual representation of trasmission clustering'
p1 <- ggpie(data = cluster_df, x = 'count', color = 'white', lab.font = "white",
            lab.pos = 'in', fill = "type") + scale_fill_npg() +
  guides(fill=guide_legend(title = "Type")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'right')
p1_pdf <- file.path(working_dir, 'clustering_visualization.pdf')
p1_png <- file.path(working_dir, 'clustering_visualization.png')
ggsave(p1_pdf, width = 7.8, height = 5.72)
ggsave(p1_png, width = 7.8, height = 5.72, dpi = 500)

####### 2nd plot: a pie chart of the number of members in each cluster

cat("[2] Second pie chart...\n")
if ('Clustered' %in% df$clustering_type) {
  new_df <- df[which(df$clustering_type == 'Clustered'),]
  s_title <- 'Visual representation of the number of members in each cluster'
  stat_df <- as.data.frame(table(new_df[,1]))
  group_df <- as.data.frame(table(stat_df[,2]))
  colnames(group_df) <- c('members_number', 'count')
  p2 <- ggpie(data = group_df, x = 'count', color = 'white', lab.font = "white",
              lab.pos = 'in', fill = "members_number") + scale_fill_jama() +
    guides(fill=guide_legend(title = "Members")) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'right')
} else {
  p2 <- ggplot() + theme_void() + geom_text(aes(0, 0, label = 'N/A')) +
    xlab(NULL) #optional, but safer in case another theme is applied later
}

p2_pdf <- file.path(working_dir, 'cluster_member_visualization.pdf')
p2_png <- file.path(working_dir, 'cluster_member_visualization.png')
ggsave(p2_pdf, width = 7.8, height = 5.72)
ggsave(p2_png, width = 7.8, height = 5.72, dpi = 500)


