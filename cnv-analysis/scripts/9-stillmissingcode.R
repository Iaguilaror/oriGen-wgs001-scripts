### CHEKPOINT
# 
# # 3. Rango de CNV por gen
# # Create the CNV matrix
# cnv_matrix_coding.df <- fix_ann_affected_genes_FR.df %>% 
#   filter( gene_biotype == "protein_coding" ) %>% 
#   select( sample, external_gene_name, CN ) %>% 
#   unique( ) %>% 
#   pivot_wider( id_cols = external_gene_name,
#                names_from = "sample",
#                values_from = "CN",
#                values_fill = 2 )
# 
# # Add columns for min, max, mean, and median
# summ_cnv_matrix_coding.df <- cnv_matrix_coding.df %>%
#   rowwise() %>%
#   mutate(
#     min_value = min(c_across(-external_gene_name)), # Exclude the gene name column
#     max_value = max(c_across(-external_gene_name)),
#     mean_value = mean(c_across(-external_gene_name)),
#     median_value = median(c_across(-external_gene_name))
#   ) %>%
#   ungroup() %>% 
#   select(external_gene_name, min_value, max_value, mean_value, median_value )
# 
# # plot histos for by gene metrics
# min_histo.p <- ggplot( data = summ_cnv_matrix_coding.df,
#                        mapping = aes( x = min_value ) ) +
#   geom_histogram( bins = 100 ) +
#   labs( title = "Dist of metrics for CNV in genes",
#         subtitle = "protein_coding" ) +
#   theme_light( base_size = 15 )
# 
# max_histo.p <- ggplot( data = summ_cnv_matrix_coding.df,
#                        mapping = aes( x = max_value ) ) +
#   geom_histogram( bins = 100 ) +
#   labs( title = "Dist of metrics for CNV in genes",
#         subtitle = "protein_coding" ) +
#   theme_light( base_size = 15 )
# 
# mean_histo.p <- ggplot( data = summ_cnv_matrix_coding.df,
#                         mapping = aes( x = mean_value ) ) +
#   geom_histogram( bins = 100 ) +
#   labs( title = "Dist of metrics for CNV in genes",
#         subtitle = "protein_coding" ) +
#   theme_light( base_size = 15 )
# 
# median_histo.p <- ggplot( data = summ_cnv_matrix_coding.df,
#                           mapping = aes( x = median_value ) ) +
#   geom_histogram( bins = 100 ) +
#   labs( title = "Dist of metrics for CNV in genes",
#         subtitle = "protein_coding" ) +
#   theme_light( base_size = 15 )
# 
# ## save the plot
# ggsave( filename = "min_histo_coding.png", plot = min_histo.p, width = 10, height = 10 )
# ## save the plot
# ggsave( filename = "max_histo_coding.png", plot = max_histo.p, width = 10, height = 10 )
# ## save the plot
# ggsave( filename = "mean_histo_coding.png", plot = mean_histo.p, width = 10, height = 10 )
# ## save the plot
# ggsave( filename = "median_histo_coding.png", plot = median_histo.p, width = 10, height = 10 )
# 
# # Create the CNV matrix
# cnv_matrix_lnc.df <- fix_ann_affected_genes_FR.df %>% 
#   filter( gene_biotype == "lncRNA" ) %>% 
#   select( sample, external_gene_name, CN ) %>% 
#   unique( ) %>% 
#   pivot_wider( id_cols = external_gene_name,
#                names_from = "sample",
#                values_from = "CN",
#                values_fill = 2 )
# 
# # Add columns for min, max, mean, and median
# summ_cnv_matrix_lnc.df <- cnv_matrix_lnc.df %>%
#   rowwise() %>%
#   mutate(
#     min_value = min(c_across(-external_gene_name)), # Exclude the gene name column
#     max_value = max(c_across(-external_gene_name)),
#     mean_value = mean(c_across(-external_gene_name)),
#     median_value = median(c_across(-external_gene_name))
#   ) %>%
#   ungroup() %>% 
#   select(external_gene_name, min_value, max_value, mean_value, median_value )
# 
# # plot histos for by gene metrics
# min_histo.p <- ggplot( data = summ_cnv_matrix_lnc.df,
#                        mapping = aes( x = min_value ) ) +
#   geom_histogram( bins = 100 ) +
#   labs( title = "Dist of metrics for CNV in genes",
#         subtitle = "lncRNA" ) +
#   theme_light( base_size = 15 )
# 
# max_histo.p <- ggplot( data = summ_cnv_matrix_lnc.df,
#                        mapping = aes( x = max_value ) ) +
#   geom_histogram( bins = 100 ) +
#   labs( title = "Dist of metrics for CNV in genes",
#         subtitle = "lncRNA" ) +
#   theme_light( base_size = 15 )
# 
# mean_histo.p <- ggplot( data = summ_cnv_matrix_lnc.df,
#                         mapping = aes( x = mean_value ) ) +
#   geom_histogram( bins = 100 ) +
#   labs( title = "Dist of metrics for CNV in genes",
#         subtitle = "lncRNA" ) +
#   theme_light( base_size = 15 )
# 
# median_histo.p <- ggplot( data = summ_cnv_matrix_lnc.df,
#                           mapping = aes( x = median_value ) ) +
#   geom_histogram( bins = 100 ) +
#   labs( title = "Dist of metrics for CNV in genes",
#         subtitle = "lncRNA" ) +
#   theme_light( base_size = 15 )
# 
# ## save the plot
# ggsave( filename = "min_histo_lnc.png", plot = min_histo.p, width = 10, height = 10 )
# ## save the plot
# ggsave( filename = "max_histo_lnc.png", plot = max_histo.p, width = 10, height = 10 )
# ## save the plot
# ggsave( filename = "mean_histo_lnc.png", plot = mean_histo.p, width = 10, height = 10 )
# ## save the plot
# ggsave( filename = "median_histo_lnc.png", plot = median_histo.p, width = 10, height = 10 )

# # Make the heatmap for CN dosage (CN*%affected gene)
# dosage_FR.df <- selected_FR.df %>% 
#   mutate( dosage = CN * feat_fracc_affected ) %>% 
#   select( sample, feat_id, dosage ) %>% 
#   left_join( x = .,
#              y = biomart.df,
#              by = c( "feat_id" = "ensembl_transcript_id_version" ) ) %>% 
#   filter(  gene_biotype == "protein_coding" | gene_biotype == "lncRNA" ) %>% 
#   mutate( external_gene_name = ifelse( is.na( external_gene_name ),
#                                        yes = ensembl_gene_id_version,
#                                        no = external_gene_name ) ) %>% 
#   unique( ) %>% 
#   arrange( sample, external_gene_name, -dosage ) %>% 
#   group_by( sample, external_gene_name ) %>% 
#   slice( 1 ) %>% 
#   ungroup( )

# make_pheatmap2.f <- function( the_data, the_biotype ){
#   
#   # prepare for pheatmap
#   wide_dosage_FR.df <- the_data %>% 
#     filter( gene_biotype == the_biotype ) %>% 
#     pivot_wider( id_cols = external_gene_name,
#                  names_from = sample,
#                  values_from = dosage, values_fill = 2  )
#   
#   # Convert the tibble to a matrix for pheatmap (excluding the external_gene_name column)
#   cnv_matrix <- as.matrix(wide_dosage_FR.df[ , -1])
#   
#   # Set row names to the external_gene_name column
#   rownames(cnv_matrix) <- wide_dosage_FR.df$external_gene_name
#   
#   # Generate the heatmap
#   pheatmap(cnv_matrix, 
#            color = colorRampPalette(c("white", "orange", "red"))(10), 
#            cluster_rows = TRUE, 
#            cluster_cols = TRUE, 
#            show_rownames = FALSE, 
#            show_colnames = FALSE,
#            main = paste( "CopyNumber dosage (CN*%affected bp)", "-", the_biotype ) )
#   
# }
# 
# png( filename = "pheatmap_coding_dosage.png", width = 1200, height = 1200, units = "px" )
# make_pheatmap2.f( the_data = dosage_FR.df, the_biotype = "protein_coding"  )
# dev.off( )
# 
# png( filename = "pheatmap_lncRNA_dosage.png", width = 1200, height = 1200, units = "px" )
# make_pheatmap2.f( the_data = dosage_FR.df, the_biotype = "lncRNA"  )
# dev.off( )

##
# vt: yo diria muestres en heatmap los que tienen mas de 3 y menos de 1000


# the_data %>% 
#   filter( gene_biotype == the_biotype ) %>% 
#   select( sample, external_gene_name, CN ) %>% 
#   unique( ) %>% 
#   pivot_wider( id_cols = external_gene_name,
#                names_from = "sample",
#                values_from = "CN",
#                values_fill = 2 )

png( filename = "pheatmap_lncRNA_100perc_affected_low-freq.png", width = 1200, height = 1200, units = "px" )
make_pheatmap.f( the_data = full_filtered_for_heatmap.df, the_biotype = "lncRNA"  )
dev.off( )
