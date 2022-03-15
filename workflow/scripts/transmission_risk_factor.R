suppressWarnings(suppressMessages(library(gtsummary)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(argparse)))
suppressWarnings(suppressMessages(library(tibble)))
suppressWarnings(suppressMessages(library(ggpubr)))

rm(list = ls())

# for table1 type of 'full'
myFun <- function(x) {
  if (x != '') {
    if (x == '<.001') {
      value <- paste0(x, '**')
    } else {
      if (as.numeric(x) < 0.01) {
        value <- paste0(x, '**')
      } else if (as.numeric(x) < 0.05) {
        value <- paste0(x, '*')
      } else {
        value <- paste0(x)
      }
    }
  } else {
    value <- ''
  }
  return(value)
}

cat(paste0('[0] Transmission risk factors investigation.\n'))
parser <- ArgumentParser(description = 'Identification of risk factors to transmission')
parser$add_argument('-k', dest = 'characteristics', required = 'TRUE',
                    help = 'Epidemiological characteristics of the cases used to identify transmission risk factors.')
parser$add_argument('-c', dest = 'clusters', required = 'TRUE',
                    help = 'Genomic clusters by whole-genome sequencing.')
parser$add_argument('-l', dest = 'list', required = 'TRUE', nargs = '+',
                    help = 'Characteristics for testing.')
parser$add_argument('-o', dest = 'output_dir', required = 'TRUE',
                    help = 'Setting output directory.')

######################## read input args ########################
args <- parser$parse_args()
char_file <- args$characteristics
clus_file <- args$clusters
c_list <- args$list
output_dir <- args$output_dir

char_df <- read.table(char_file, sep = '\t', header = T, check.names = F, comment.char = "#", stringsAsFactors = F)
clus_df <- read.table(clus_file, sep = '\t', header = T, check.names = F, comment.char = "#", stringsAsFactors = F)
setwd(output_dir)
total <- merge(clus_df, char_df, by = 'sample')
###############################################################################################
cat(paste0('[1] Epidemiological characteristics: [', paste(c_list, collapse = ', '), '].\n'))
valid_cols <- c('clustering')
for (i in c_list) {
  if (i %in% colnames(total)) {
    if (class(total[,i]) == 'character') total[,i] <- as.factor(total[,i])
    valid_cols <- c(valid_cols, i)
  } else {
    print(paste0('Warning, [', i, '] is not in the characters file.'))
  }
}
if (length(valid_cols) == 0) {
  stop("Found no character in the input file.")
}
#
cat(paste0('[2] Using univariable test to compute associated significance.'))

total <- total %>% add_column(clustering = ifelse(.$clustering_type == 'Clustered', 1, 0))
select_data <- total %>% select(all_of(valid_cols))
write.table(x = select_data, file = 'characteristic_data.txt', sep = '\t', quote = F, fileEncoding = 'UTF-8', row.names = F)

################# perform univariate regression test #######################
set_gtsummary_theme(theme_gtsummary_journal("jama"))

tbl <- select_data %>%
  tbl_uvregression(
    method = glm,
    y = 'clustering',
    method.args = list(family = binomial),
    exponentiate = T,
    hide_n = T,
    pvalue_fun = ~style_pvalue(.x, digits = 2)
  ) %>%
  add_global_p(quiet = T) %>%  # add global p-value
  bold_p() %>%        # bold p-values under a given threshold (default 0.05)
  bold_p(t = 0.05) %>% # now bold q-values under the threshold of 0.10
  bold_labels() %>%
  modify_table_body(
    ~ .x %>%
      dplyr::mutate(u_event = n_obs - n_event) %>%
      dplyr::relocate(c(n_event, u_event), .before = estimate)
  ) %>%
  # assigning header labels
  modify_header(n_event = "**Clustered**", u_event = "**Unique**") %>%
  modify_fmt_fun(c(n_event, u_event) ~ style_number)
gt::gtsave(as_gt(tbl), file = 'univariable_regression_result.html')
gt::gtsave(as_gt(tbl), file = 'univariable_regression_result.rtf')
