# Load required packages
pacman::p_load("vroom", "dplyr", "tidyr")

# Read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# Debugging mode (uncomment for testing)
# args[1] <- "gt_matrix.tmp"
# args[2] <- "random100000_variants_clean_vep_loftee_only.tsv"

# Assign arguments
vartsv <- args[1]
output_file <- args[2]

# Load data efficiently
onlycounts.df <- vroom(file = vartsv, show_col_types = FALSE) %>%
  as_tibble()

# Step 1: Keep only HC LoF variants
#hc_only.df <- onlycounts.df %>%
#  filter(LoF == "HC")

# Step 2: Recode genotypes
recode_genotype <- function(gt) {
  case_when(
    gt == "0/0" ~ 0,
    gt %in% c("0/1", "1/0") ~ 1,
    gt == "1/1" ~ 2,
    TRUE ~ NA_real_
  )
}

# Apply recoding to genotype columns only (dynamically detected)
sample_cols <- which(!names(onlycounts.df) %in% c("TYPE", "LoF"))  # Assumes TYPE & LoF are the only non-genotype cols

recoded.df <- onlycounts.df %>%
  mutate(across(all_of(sample_cols), recode_genotype))

# Step 3: Count genotype categories per sample and TYPE
genotype_summary.df_1 <- recoded.df %>%
  pivot_longer(cols = -TYPE:-LoF, names_to = "sample", values_to = "gt") %>%
  count(TYPE, sample, gt) %>%
  pivot_wider(
    names_from = gt,
    values_from = n,
    values_fill = 0,
    names_sort = TRUE
  ) 

# check if col2 does NOT exist
if( ! "2" %in% colnames( genotype_summary.df_1 ) ) {
  genotype_summary.df_1 <- genotype_summary.df_1 %>% 
  mutate( "2" = 0 )
}

genotype_summary.df <-  genotype_summary.df_1 %>%
  rename(
    hom_ref = `0`,
    het = `1`,
    hom_LoF = `2`
  )

# Export result
vroom_write(genotype_summary.df, file = output_file, delim = "\t", col_names = TRUE)
