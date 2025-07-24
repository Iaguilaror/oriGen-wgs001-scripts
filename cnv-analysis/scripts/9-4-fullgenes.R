# load packages
pacman::p_load( "vroom", "dplyr", "tidyr", "ggplot2", "stringr", "ggsci", "tibble",
                "ggExtra", "patchwork", "scales", "pheatmap", "ggvenn" )

# load functions for project
source( file = "9-functions.R" )

# load data
load( file = "forfullgenes.RData" )

full_fix_ann_affected_genes_FR.df <- ann_affected_genes_FR.df %>% 
  select( sample, CN, type,
          ensembl_gene_id_version, external_gene_name, gene_biotype, feat_fracc_affected ) %>% 
  mutate( external_gene_name = ifelse( is.na( external_gene_name ),
                                       yes = ensembl_gene_id_version,
                                       no = external_gene_name ) ) %>% 
  filter( feat_fracc_affected == 1  )

full_base_metrics.df <- full_fix_ann_affected_genes_FR.df %>% 
  select( CN, sample, external_gene_name, gene_biotype  ) %>% 
  group_by( external_gene_name, gene_biotype ) %>% 
  summarise( min_CN = min( CN ),
             max_CN = max( CN ),
             mean_CN = mean( CN ),
             meadian_CN = median( CN ),
             samples = n( ) ) %>% 
  mutate( delta = max_CN - min_CN ) %>% 
  ungroup( )

# get number of fully dels and dups by gene
full_genes_by_cnv_type.df <- full_fix_ann_affected_genes_FR.df %>% 
  group_by( external_gene_name, type ) %>% 
  summarise( n_samples = n( ) )

# write the sumary
write.table( x = full_genes_by_cnv_type.df,
             file = "full_genes_by_cnv_type.tsv",
             append = F, quote = F, sep = "\t", row.names = F, col.names = T )

# Work with genes in more than 15 samples
genes_morethen15samples.df <- full_base_metrics.df %>% 
  filter( samples > 15 )

range_coding.p <- do_rangeplot.f( the_data = genes_morethen15samples.df,
                                  the_biotype = "protein_coding" )

range_lncrna.p <- do_rangeplot.f( the_data = genes_morethen15samples.df,
                                  the_biotype = "lncRNA" )

## save the plot
ggsave( filename = "range_coding.png", plot = range_coding.p, width = 9, height = 9 )
## save the plot
ggsave( filename = "range_lncrna.png", plot = range_lncrna.p, width = 9, height = 9 )

# save the data for those plots
write.table( x = genes_morethen15samples.df,
             file = "100percent_affected_genes_inmorethan15samples.tsv",
             append = F, quote = F, sep = "\t",
             row.names = F, col.names = T )

# summarise
full_summ_gene_by_sample.df <- full_fix_ann_affected_genes_FR.df %>% 
  group_by( sample, type, gene_biotype ) %>% 
  summarise( n = n( ) ) %>% 
  ungroup( )

# Function to create scatterplot with marginal boxplots
ncvns_fullaffected_genes_scatter_coding.p <- full_summ_gene_by_sample.df %>% filter( gene_biotype == "protein_coding") %>% 
  pivot_wider( id_cols = sample,
               names_from = type,
               values_from =  n, values_fill = 0 ) %>% 
  create_scatter_with_marginals.f( data = ., x_col = "DEL", y_col = "DUP",
                                   the_fill = "darkgreen", the_shape = 20,  the_alpha = 0.2,
                                   the_size = 3, base_size = 15,
                                   the_title = "Number Genes Affected by CNV",
                                   the_subtitle = "protein_coding",
                                   x_title = "n DELetions",
                                   y_title = "n DUPlications",
                                   the_caption = "For genes 100% DELeted or DUPlicated",
                                   the_hline = NA )

ncvns_fullaffected_genes_scatter_lncrna.p <- full_summ_gene_by_sample.df %>% filter( gene_biotype == "lncRNA") %>% 
  pivot_wider( id_cols = sample,
               names_from = type,
               values_from =  n, values_fill = 0 ) %>% 
  create_scatter_with_marginals.f( data = ., x_col = "DEL", y_col = "DUP",
                                   the_fill = "brown", the_shape = 20,  the_alpha = 0.2,
                                   the_size = 3, base_size = 15,
                                   the_title = "Number Genes Affected by CNV",
                                   the_subtitle = "lncRNA",
                                   x_title = "n DELetions",
                                   y_title = "n DUPlications",
                                   the_caption = "For genes 100% DELeted or DUPlicated",
                                   the_hline = NA )

# save plots
ggsave( filename = "ncvns_fullaffected_genes_scatter_coding.png",
        plot = ncvns_fullaffected_genes_scatter_coding.p, width = 7, height = 7 )

# save plots
ggsave( filename = "ncvns_fullaffected_genes_scatter_lncrna.png",
        plot = ncvns_fullaffected_genes_scatter_lncrna.p, width = 7, height = 7 )

# summarise for report
summ_ngenes.f( the_data = full_summ_gene_by_sample.df, the_biotype = "protein_coding", the_cnvtype = "DUP" )
summ_ngenes.f( the_data = full_summ_gene_by_sample.df, the_biotype = "protein_coding", the_cnvtype = "DEL" )

summ_ngenes.f( the_data = full_summ_gene_by_sample.df, the_biotype = "lncRNA", the_cnvtype = "DUP" )
summ_ngenes.f( the_data = full_summ_gene_by_sample.df, the_biotype = "lncRNA", the_cnvtype = "DEL" )

# Make a plot to compare sets
fullgenes_barvenn_coding.list <- genes_morethen15samples.df %>% 
  filter( gene_biotype == "protein_coding" ) %>% 
  barvenn.f( the_data = .,
             the_subtitle = "protein_coding - genes in more than 1% of samples" )

fullgenes_barvenn_coding.p <- fullgenes_barvenn_coding.list$the_plot

fullgenes_barvenn_coding.df <- fullgenes_barvenn_coding.list$the_df

write.table( x = fullgenes_barvenn_coding.df, file = "fullgenes_barvenn_coding.tsv",
             append = F, quote = F, sep = "\t", row.names = F, col.names = T )

fullgenes_barvenn_lncrna.list <- genes_morethen15samples.df %>% 
  filter( gene_biotype == "lncRNA" ) %>% 
  barvenn.f( the_data = .,
             the_subtitle = "lncRNA - genes in more than 1% of samples" )

fullgenes_barvenn_lncrna.p <- fullgenes_barvenn_lncrna.list$the_plot

fullgenes_barvenn_lncrna.df <- fullgenes_barvenn_lncrna.list$the_df

write.table( x = fullgenes_barvenn_lncrna.df, file = "fullgenes_barvenn_lncrna.tsv",
             append = F, quote = F, sep = "\t", row.names = F, col.names = T )

# save plots
ggsave( filename = "fullgenes_barvenn_coding.png", bg = "white",
        plot = fullgenes_barvenn_coding.p, width = 10, height = 7 )

# save plots
ggsave( filename = "fullgenes_barvenn_lncrna.png", bg = "white",
        plot = fullgenes_barvenn_lncrna.p, width = 10, height = 7 )

# Do heatmaps
# get max number of samples
max_samples.v <- ann_affected_genes_FR.df %>% 
  pull( sample ) %>% 
  unique( ) %>% 
  length( )

# get filtered external_gene_name
full_filtered_for_heatmap.v <- full_base_metrics.df %>% 
  filter( samples >= max_samples.v * 0.01   ) %>% 
  filter( samples <= max_samples.v * 0.05   ) %>% 
  pull( external_gene_name ) %>% 
  unique( )

full_filtered_for_heatmap.df <- full_fix_ann_affected_genes_FR.df %>% 
  filter( external_gene_name %in% full_filtered_for_heatmap.v )

dev.off( )

png( filename = "pheatmap_protein_100perc_affected_low-freq.png", width = 1200, height = 1000, units = "px" )
make_pheatmap.f( the_data = full_filtered_for_heatmap.df, the_biotype = "protein_coding"  )
dev.off( )

png( filename = "pheatmap_ncrna_100perc_affected_low-freq.png", width = 1200, height = 1000, units = "px" )
make_pheatmap.f( the_data = full_filtered_for_heatmap.df, the_biotype = "lncRNA"  )
dev.off( )

# Now high freqs
# get filtered external_gene_name
full_filtered_for_heatmap.v <- full_base_metrics.df %>% 
  filter( samples >= max_samples.v * 0.05   ) %>% 
  pull( external_gene_name ) %>% 
  unique( )

full_filtered_for_heatmap.df <- full_fix_ann_affected_genes_FR.df %>% 
  filter( external_gene_name %in% full_filtered_for_heatmap.v )

# save the plots
png( filename = "pheatmap_protein_100perc_affected_common-freq.png", width = 1200, height = 1000, units = "px" )
make_pheatmap.f( the_data = full_filtered_for_heatmap.df, the_biotype = "protein_coding"  )
dev.off( )

png( filename = "pheatmap_ncrna_100perc_affected_common-freq.png", width = 1200, height = 1000, units = "px" )
make_pheatmap.f( the_data = full_filtered_for_heatmap.df, the_biotype = "lncRNA"  )
dev.off( )


# POSANALYSIS for easy habdling of results for paper #####
fullgenes_barvenn_coding.df.annotated <- left_join( x = fullgenes_barvenn_coding.df,
                                                    y = genes_morethen15samples.df,
                                                    by = c("genes" = "external_gene_name" ) )

write.table( x = fullgenes_barvenn_coding.df.annotated,
             file = "fullgenes_barvenn_coding.df.annotated.tsv",
             append = F, quote = F, sep = "\t",
             row.names = F, col.names = T )
