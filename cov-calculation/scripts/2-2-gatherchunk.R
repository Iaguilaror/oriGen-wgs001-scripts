#!/usr/bin/env Rscript

# -------------------------------
# Load required packages
# -------------------------------
#if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyr, dplyr, vroom, purrr)

# -------------------------------
# List input files matching pattern
# -------------------------------
input_files <- list.files(
  path = ".",
  pattern = "_samtools_mean_dp_all_cov_for_cnv.tsv$",
  full.names = TRUE
)

# -------------------------------
# Read and join all files by key columns
# -------------------------------
joined_df <- input_files %>%
  set_names() %>%  # use file names as names for each element
  map(~ vroom(.x, show_col_types = FALSE)) %>%
  reduce(full_join, by = c("chromosome", "start", "end", "gene_name", "exon_rank", "transcript_id"))

# -------------------------------
# Write the combined result to a compressed file
# -------------------------------
output_file <- "gathered_mean_dp_all_cov_for_cnv.tsv.gz"
vroom_write(joined_df, output_file, delim = "\t", na = "")
