# load packages
pacman::p_load( "vroom", "dplyr", "tidyr", "ggplot2", "stringr", "ggsci",
                "ggExtra", "patchwork", "scales", "pheatmap", "ggvenn" )

# load functions
source( file = "9-functions.R" )

# load FR data
all_FR.df <- vroom( file = "allsamples.FR.bed", col_names = FALSE ) %>% 
  as_tibble( )

selected_FR.df <- all_FR.df %>% 
  select( X10, X9, X4, X8, X11, X12 ) %>% 
  rename( sample = 1,
          cnv_id = 2,
          feat_id = 3,
          CN = 4,
          feat_bp_affected = 5,
          feat_fracc_affected = 6 ) %>% 
  mutate(type = case_when(
    str_detect(cnv_id, "del") ~ "DEL",
    str_detect(cnv_id, "dup") ~ "DUP",
    TRUE ~ "other"  ) )

# get uniq ids by type
by_sample_n_cnv.df <- selected_FR.df %>% 
  select( sample, cnv_id, type ) %>% 
  unique( ) %>% 
  group_by( sample, type ) %>% 
  summarise( n = n( ) ) %>% 
  ungroup( )

# plot DEL
delbysample.p <- create_bars_ordered_bysample.f( the_data = by_sample_n_cnv.df, col_to_analyze = n,
                                                 the_type = "DEL", the_fill = "tomato",
                                                 the_title =  "CNV DELetion distribution in oriGEN phase 1",
                                                 the_x_title = "Sample", the_y_title = "Number of DELs" )

dupbysample.p <- create_bars_ordered_bysample.f( the_data = by_sample_n_cnv.df, col_to_analyze = n,
                                                 the_type = "DUP", the_fill = "steelblue",
                                                 the_title =  "CNV DUPlication distribution in oriGEN phase 1",
                                                 the_x_title = "Sample", the_y_title = "Number of DUPs" )

# save the by sample plots
ggsave( filename = "n_delbysample.png", plot = delbysample.p, width = 10, height = 7 )
ggsave( filename = "n_dupbysample.png", plot = dupbysample.p, width = 10, height = 7 )

# get summary for numbers
by_sample_n_cnv.df %>% 
  filter( type == "DEL" ) %>% 
  summary( )

by_sample_n_cnv.df %>% 
  filter( type == "DUP" ) %>% 
  summary( )

# create a scatter with number of dels ####
# scatter to compare bp DEL vs bp DUP
# Reshape the dataframe to wide format with separate columns for DEL and DUP
scatter_data1 <- by_sample_n_cnv.df %>%
  pivot_wider( names_from = type, values_from = n )

# Create the scatterplot
n_cnv_scatterside.p <- create_scatter_with_marginals.f( data = scatter_data1,
                                                        x_col = "DEL", y_col = "DUP",
                                                        the_fill = "limegreen", the_shape = 21,
                                                        the_size = 3, base_size = 15,
                                                        the_title = "Distribution of Number of CNVs",
                                                        the_subtitle = "each dot is an oriGEN sample",
                                                        x_title = "Number of DEL", y_title = "Number of DUP" )

## save the plot
ggsave( filename = "n_cnv_scatterside.png", plot = n_cnv_scatterside.p, width = 10, height = 10 )

## Summarise the total number of bps affected by sample by type ####
by_sample_totalbp_cnv.df <- selected_FR.df %>%
  select( sample, feat_id, feat_bp_affected, type ) %>%
  unique( ) %>%
  group_by( sample, type ) %>%
  summarise( total_affected_bp = sum( feat_bp_affected )  ) %>%
  ungroup( )

bp_delbysample.p <- create_bars_ordered_bysample.f( the_data = by_sample_totalbp_cnv.df, col_to_analyze = total_affected_bp,
                                                    the_type = "DEL", the_fill = "tomato", the_alpha = 0.5,
                                                    the_title =  "Total Affected Base Pairs - DELetions",
                                                    the_x_title = "Sample", the_y_title = "Total DELETED bp" )

bp_dupbysample.p <- create_bars_ordered_bysample.f( the_data = by_sample_totalbp_cnv.df, col_to_analyze = total_affected_bp,
                                                    the_type = "DUP", the_fill = "steelblue", the_alpha = 0.5,
                                                    the_title =  "Total Affected Base Pairs - DUPlications",
                                                    the_x_title = "Sample", the_y_title = "Total DUPLICATED bp" )

# save the by sample plots
ggsave( filename = "bp_delbysample.png", plot = bp_delbysample.p, width = 10, height = 7 )
ggsave( filename = "bp_dupbysample.png", plot = bp_dupbysample.p, width = 10, height = 7 )

# get summary for numbers
by_sample_totalbp_cnv.df %>% 
  filter( type == "DEL" ) %>% 
  summary( )

by_sample_totalbp_cnv.df %>% 
  filter( type == "DUP" ) %>% 
  summary( )

# scatter to compare bp DEL vs bp DUP
# Reshape the dataframe to wide format with separate columns for DEL and DUP
scatter_data <- by_sample_totalbp_cnv.df %>%
  pivot_wider(names_from = type, values_from = total_affected_bp)

bp_scatter.p <- create_scatter_with_marginals.f( data = scatter_data,
                                                 x_col = "DEL", y_col = "DUP",
                                                 the_fill = "orange", the_shape = 21, the_alpha = 0.5,
                                                 the_size = 3, base_size = 15,
                                                 the_title = "Distribution of Total Affected BP by CNVs",
                                                 the_subtitle = "each dot is an oriGEN sample",
                                                 x_title = "log10 Total Affected BP (DEL)",
                                                 y_title = "log10 Total Affected BP (DUP)" )

## save the plot
ggsave( filename = "bp_scatter.png", plot = bp_scatter.p, width = 14, height = 12 )

# Let's Answer the following
# 1. N de variantes CNV x muestra
# it's in n_scatter_side.png
# Save the dataframe: scatter_data1
write.table( x = scatter_data1 , file = "n_CNV_by_sample.tsv",
             append = F, quote = F,
             sep = "\t", row.names = F, col.names = T )

#### POSSIBLE CHECKPOINT ####
