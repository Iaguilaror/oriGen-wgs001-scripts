# load pkgs
pacman::p_load( "vroom", "dplyr", "tidyr", "ggplot2" )

# Find all .tsv files recursively in the current directory
tsv_files <- list.files( path = ".",
                         pattern = "\\.tsv$",
                         full.names = TRUE,
                         recursive = TRUE )

# Read and combine all .tsv files into a single dataframe
combined_df <- vroom( tsv_files, id = "source_file", show_col_types = F ) %>% 
  as_tibble( )

str( combined_df )

summary_by_sample_type <- combined_df %>%
  group_by(sample, TYPE) %>%
  summarise(
    hom_ref = sum(hom_ref, na.rm = TRUE),
    het = sum(het, na.rm = TRUE),
    hom_LoF = sum(hom_LoF, na.rm = TRUE),
    .groups = "drop"
  )

str( summary_by_sample_type )

# Reshape to long format for easier plotting
long_summary <- summary_by_sample_type %>%
  pivot_longer(
    cols = c(hom_ref, het, hom_LoF),
    names_to = "genotype",
    values_to = "count"
  )

# Plot boxplot
boxplots.p <- long_summary %>%
	filter(genotype != "hom_ref" ) %>%
ggplot( data = ., mapping = aes(x = genotype, y = count, fill = genotype)) +
  geom_boxplot() +
  facet_wrap(~TYPE, ncol = 1 ) +
  theme_minimal( base_size = 15 ) +
  labs(
    title = "Genotype counts per sample by variant TYPE",
    x = "Genotype category",
    y = "Count per sample"
  )

ggsave( filename = "n_HC_LoF_variants.svg",
        plot = boxplots.p,
        width = 10, height = 7 )

# resume data
summary_stats <- summary_by_sample_type %>%
  pivot_longer(
    cols = c(hom_ref, het, hom_LoF),
    names_to = "genotype",
    values_to = "count"
  ) %>%
  group_by(TYPE, genotype) %>%
  summarise(
    mean   = mean(count, na.rm = TRUE),
    median = median(count, na.rm = TRUE),
    sd     = sd(count, na.rm = TRUE),
    se     = sd / sqrt(n()),
    min    = min(count, na.rm = TRUE),
    max    = max(count, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  )

# save the data
vroom_write( x = summary_stats, file = "summary_stats_LoF_allsamples.tsv", "\t" )
