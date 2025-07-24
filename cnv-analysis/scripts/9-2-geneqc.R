# load packages
pacman::p_load( "vroom", "dplyr", "tidyr", "ggplot2", "stringr", "ggsci",
                "ggExtra", "patchwork", "scales", "pheatmap", "ggvenn" )

# load functions of project
source( file = "9-functions.R" )

# define params
gene_distance_to_gap <- 30000
gene_length_cut <- 0.5

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

# write data for figures
write.table( x = selected_FR.df,
             file = "selected_FR.tsv",
             append = F, quote = F,
             sep = "\t",
             row.names = F, col.names = T )

##### Start by gene #####
cn_by_fracc_del.p <- create_scatter_with_marginals.f( data = filter( selected_FR.df, type == "DEL" ),
                                                      x_col = "CN", y_col = "feat_fracc_affected", the_fill = "tomato",
                                                      the_shape = 21, the_alpha = 0.5,
                                                      the_size = 3, base_size = 15,
                                                      the_title = "Distribution of DELetions", the_subtitle = "",
                                                      x_title = "Copy Number (CN)",
                                                      y_title = "Length of Gene Affected (feat_fracc_affected)",
                                                      the_hline = gene_length_cut )

cn_by_fracc_dup.p <- create_scatter_with_marginals.f( data = filter( selected_FR.df, type == "DUP" ),
                                                      x_col = "CN", y_col = "feat_fracc_affected", the_fill = "skyblue",
                                                      the_shape = 21, the_alpha = 0.5,
                                                      the_size = 3, base_size = 15,
                                                      the_title = "Distribution of DUPlications", the_subtitle = "",
                                                      x_title = "Copy Number (CN)",
                                                      y_title = "Length of Gene Affected (feat_fracc_affected)",
                                                      the_hline = gene_length_cut  )

## save the plot
ggsave( filename = "cn_by_fracc_del.png", plot = cn_by_fracc_del.p, width = 10, height = 10 )
## save the plot
ggsave( filename = "cn_by_fracc_dup.png", plot = cn_by_fracc_dup.p, width = 10, height = 10 )

# get numbers of DEL
selected_FR.df %>% 
  filter( type == "DEL" ) %>% 
  nrow( )

# get numbers of DEL
selected_FR.df %>% 
  filter( type == "DUP" ) %>% 
  nrow( )

# now filter, only genes(features) affected by at least 50% of its length
affected_genes_FR.df <- selected_FR.df %>% 
  filter( feat_fracc_affected >= gene_length_cut )

# get numbers of DEL
affected_genes_FR.df %>% 
  filter( type == "DEL" ) %>% 
  nrow( )

# get numbers of DEL
affected_genes_FR.df %>% 
  filter( type == "DUP" ) %>% 
  nrow( )

# get gene distance to telomere and centromere
# read telomere
telomere_base.df <- vroom( file = "telomere_positions.txt", col_names = FALSE ) %>% 
  as_tibble( ) %>% 
  select( 1, 3, 4, 6 ) %>% 
  rename( chrom = 1,
          chromStart = 2,
          chromEnd = 3,
          gaptype = 4 ) %>% 
  group_by( chrom ) %>% 
  mutate( gaptype = paste0( gaptype, 1:2 ) ) %>% 
  ungroup( ) %>% 
  pivot_longer( cols = c( chromStart, chromEnd ),
                names_to =  "vars",
                values_to = "pos" ) %>% 
  mutate( vars = str_remove( string = vars, patter = "chrom" ) ) %>% 
  mutate( newvar = paste( gaptype, vars, sep = "-" ) ) %>% 
  select( -gaptype, -vars )

# read centromere
centromere.df <- vroom( file = "centromeres.tsv" ) %>% 
  as_tibble( ) %>% 
  mutate( chrom = str_remove( string = chrom, pattern = "chr" ) ) %>% 
  select( chrom, chromStart, chromEnd,  gieStain ) %>% 
  rename( gaptype = gieStain ) %>% 
  mutate( gaptype = "centromere" ) %>% 
  group_by( chrom ) %>% 
  mutate( gaptype = paste0( gaptype, 1:2 ) ) %>% 
  ungroup( ) %>% 
  pivot_longer( cols = c( chromStart, chromEnd ),
                names_to =  "vars",
                values_to = "pos" ) %>% 
  mutate( vars = str_remove( string = vars, patter = "chrom" ) ) %>% 
  mutate( newvar = paste( gaptype, vars, sep = "-" ) ) %>% 
  select( -gaptype, -vars )

# join cents and teloms
all_gaps.df <- bind_rows( telomere_base.df, centromere.df ) %>% 
  mutate( chrom = factor( chrom, levels = c( 1:22, "X", "Y" ) ) ) %>% 
  arrange( chrom, pos ) %>% 
  pivot_wider( id_cols = chrom, names_from = "newvar", values_from = "pos" )

all_gaps_size.df <- all_gaps.df %>% 
  group_by( chrom ) %>% 
  mutate( tel_size = max( c(`telomere1-End` - `telomere1-Start`), c(`telomere2-End` - `telomere2-Start`)  ) ) %>%
  mutate( cen_size = max( `centromere2-End` - `centromere1-Start`) ) %>% 
  ungroup( ) %>% 
  select( chrom, tel_size, cen_size )

# lets annotate the genes
biomart.df <- vroom( file = "full_biomaRt_used.tsv" )

# annotate genes
ann_affected_genes_FR.df <- affected_genes_FR.df %>% 
  left_join( x = .,
             y = biomart.df,
             by = c( "feat_id" = "ensembl_transcript_id_version" ) ) %>% 
  filter(  gene_biotype == "protein_coding" | gene_biotype == "lncRNA" )

## Here filter by distance to gaps (teloms and centro)
affected_genes_simplified.df <- ann_affected_genes_FR.df %>% 
  select( feat_id, chromosome_name, start_position, end_position ) %>% 
  unique( ) %>% 
  rename( chrom = chromosome_name,
          chromStart = start_position,
          chromEnd = end_position )

#####
# Step 1: Add a gap name column and reshape the gaps
gap_regions <- all_gaps.df %>%
  pivot_longer(
    cols = starts_with("telomere") | starts_with("centromere"),
    names_to = "gap_name",
    values_to = "position"
  ) %>%
  separate(gap_name, into = c("gap_name", "boundary"), sep = "-") %>%
  pivot_wider(names_from = boundary, values_from = position) %>%
  rename(gap_start = Start, gap_end = End)

# Step 2: Calculate distances to all gaps
annotated_genes <- affected_genes_simplified.df %>%
  rowwise() %>%
  mutate(
    closest_gap_info = list(
      gap_regions %>%
        filter(chrom == chrom) %>%
        mutate(
          distance = pmin(abs(chromStart - gap_end), abs(chromEnd - gap_start),
                          abs(chromEnd - gap_end), abs(chromStart - gap_start))
        ) %>%
        arrange(distance) %>%
        slice(1) %>%
        select(gap_name, distance)
    ),
    closest_gap_name = closest_gap_info[[1, "gap_name"]],
    closest_distance = closest_gap_info[[1, "distance"]]
  ) %>%
  select(-closest_gap_info) %>%
  ungroup()

# plot closests distances
dotplot.p <- ggplot( data = annotated_genes ) +
  geom_point( mapping = aes( x = closest_distance,
                             y = chrom,
                             fill = closest_gap_name ),
              shape = 21, alpha = 0.3 ) +
  geom_point( data = all_gaps_size.df,
              mapping = aes( x = cen_size, y = chrom, color = "Centromere Length" ),
              shape = 13, size = 3, alpha = 0.5 ) +
  scale_x_continuous( trans = "log10",
                      breaks = scales::breaks_log(10),
                      labels = scales::label_number(accuracy = 1) ) +
  scale_fill_futurama( ) +
  # guides( color = "none" ) + # Remove the color legend
  labs( title = "CNV affected Genes",
        x = "Closest distance to Gap (bp)",
        y = "Chromosome",
        fill = "Nearest Gap",
        color = "" ) +
  theme_classic( ) +
  theme( legend.position = "top" )

distance_histo.p <- ggplot( data = annotated_genes,
                            mapping = aes( x = closest_distance ) ) +
  geom_vline( xintercept = gene_distance_to_gap, alpha = 0.5, lty = "dashed", color = "tomato" ) +
  geom_histogram( alpha = 0.5, bins = 100 ) +
  scale_x_continuous( trans = "log10",
                      breaks = scales::breaks_log(10),
                      labels = scales::label_number(accuracy = 1) ) +
  theme_minimal( ) +
  theme( axis.title.x = element_blank( ),
         panel.grid = element_blank( ) )

# Combine the plots in a column
distance_panel.p <- dotplot.p / distance_histo.p + 
  plot_layout(heights = c(0.9, 0.1)) # Set relative heights: 90% for dotplot.p, 10% for distance_histo.p

ggsave( filename = "distance_panel.png", plot = distance_panel.p, width = 12, height = 7 )

#####
# CUT genes by distance to closests gap
filtered_annotated_genes.v <- annotated_genes %>% 
  filter( closest_distance >= gene_distance_to_gap ) %>% 
  pull( feat_id )

ann_affected_genes_FR.df <- ann_affected_genes_FR.df %>% 
  filter( feat_id %in% filtered_annotated_genes.v )
##### DONE
######

# now... number of genes affected by CN by CN type x biotype
fix_ann_affected_genes_FR.df <- ann_affected_genes_FR.df %>% 
  select( sample, CN, type,
          ensembl_gene_id_version, external_gene_name, gene_biotype ) %>% 
  mutate( external_gene_name = ifelse( is.na( external_gene_name ),
                                       yes = ensembl_gene_id_version,
                                       no = external_gene_name ) )

uniq_fix_ann_affected_genes_FR.df <- fix_ann_affected_genes_FR.df %>% 
  select( sample, type, external_gene_name, gene_biotype ) %>% 
  unique( )

# summarize by sample
summ_gene_by_sample.df <- uniq_fix_ann_affected_genes_FR.df %>% 
  group_by( sample, type, gene_biotype ) %>% 
  summarise( n = n( ) ) %>% 
  ungroup( )

# create plots
ngenes_affected_coding.p <- summ_gene_by_sample.df %>% filter( gene_biotype == "protein_coding" ) %>% 
  pivot_wider( id_cols = sample,
               names_from = type,
               values_from =  n, values_fill = 0 ) %>% 
  create_scatter_with_marginals.f( data = .,
                                   x_col = "DEL", y_col = "DUP",
                                   the_fill = "darkseagreen",
                                   the_shape = 21, the_alpha = 0.5,
                                   the_size = 3, base_size = 15,
                                   the_title =  "Number of Genes Affected by CNVs",
                                   the_subtitle = "protein_coding",
                                   x_title = "DEL",
                                   y_title = "DUP",
                                   the_caption = paste( "only considering genes affected in at least", percent( gene_length_cut ), "of its length" ) )

ngenes_affected_ncrna.p <- summ_gene_by_sample.df %>% filter( gene_biotype == "lncRNA" ) %>% 
  pivot_wider( id_cols = sample,
               names_from = type,
               values_from =  n, values_fill = 0 ) %>% 
  create_scatter_with_marginals.f( data = .,
                                   x_col = "DEL", y_col = "DUP",
                                   the_fill = "brown",
                                   the_shape = 21, the_alpha = 0.5,
                                   the_size = 3, base_size = 15,
                                   the_title =  "Number of Genes Affected by CNVs",
                                   the_subtitle = "lncRNA",
                                   x_title = "DEL",
                                   y_title = "DUP",
                                   the_caption = paste( "only considering genes affected in at least", percent( gene_length_cut ), "of its length" ) )

## save the plot
ggsave( filename = "ngenes_affected_coding.png", plot = ngenes_affected_coding.p, width = 10, height = 10 )
## save the plot
ggsave( filename = "ngenes_affected_ncrna.png", plot = ngenes_affected_ncrna.p, width = 10, height = 10 )

summ_ngenes.f( the_data = summ_gene_by_sample.df,
               the_biotype = "protein_coding",
               the_cnvtype = "DEL" )

summ_ngenes.f( the_data = summ_gene_by_sample.df,
               the_biotype = "protein_coding",
               the_cnvtype = "DUP" )

summ_ngenes.f( the_data = summ_gene_by_sample.df,
               the_biotype = "lncRNA",
               the_cnvtype = "DEL" )

summ_ngenes.f( the_data = summ_gene_by_sample.df,
               the_biotype = "lncRNA",
               the_cnvtype = "DUP" )

# plot number of protein coding genes affected by CNV
create_col_gene_by_sample.f <- function( the_data, the_biotype, the_cnv_type  ) {
  
  the_data %>% 
    filter( gene_biotype == the_biotype ) %>% 
    filter( type == the_cnv_type ) %>% 
    arrange( n ) %>% 
    mutate( sample = factor( sample, levels = .$sample ) ) %>% 
    ggplot( data = ., 
            mapping = aes( x = sample, y = n ) ) +
    geom_col( ) +
    labs( title = "Number of genes affected by CNV",
          subtitle = paste( the_cnv_type, the_biotype, sep = " - " ),
          caption = paste( "only considering genes affected in at least", percent( gene_length_cut ), "of its length" ), 
          y = "number of genes",
    ) +
    theme_classic( base_size = 15 ) +
    theme( axis.text.x = element_blank( ),
           axis.ticks.x = element_blank( ) )
}

ngenes_coding_del.p <- create_col_gene_by_sample.f( the_data = summ_gene_by_sample.df,
                                                    the_biotype = "protein_coding",
                                                    the_cnv_type = "DEL" )

ngenes_coding_dup.p <- create_col_gene_by_sample.f( the_data = summ_gene_by_sample.df,
                                                    the_biotype = "protein_coding",
                                                    the_cnv_type = "DUP" )

ngenes_lnc_del.p <- create_col_gene_by_sample.f( the_data = summ_gene_by_sample.df,
                                                 the_biotype = "lncRNA",
                                                 the_cnv_type = "DEL" )

ngenes_lnc_dup.p <- create_col_gene_by_sample.f( the_data = summ_gene_by_sample.df,
                                                 the_biotype = "lncRNA",
                                                 the_cnv_type = "DUP" )


## save the plot
ggsave( filename = "ngenes_coding_del.png", plot = ngenes_coding_del.p, width = 10, height = 10 )
## save the plot
ggsave( filename = "ngenes_coding_dup.png", plot = ngenes_coding_dup.p, width = 10, height = 10 )
## save the plot
ggsave( filename = "ngenes_lnc_del.png", plot = ngenes_lnc_del.p, width = 10, height = 10 )
## save the plot
ggsave( filename = "ngenes_lnc_dup.png", plot = ngenes_lnc_dup.p, width = 10, height = 10 )

# 2. Genes que más/menos varían en CNVs
summ_by_gene.df <- uniq_fix_ann_affected_genes_FR.df %>% 
  group_by( external_gene_name, type, gene_biotype ) %>% 
  summarise( n = n( ) )

## add % of affection
summ_aff_by_gene.df <- ann_affected_genes_FR.df %>% 
  mutate( external_gene_name = ifelse( is.na( external_gene_name ),
                                       yes = ensembl_gene_id_version,
                                       no = external_gene_name ) ) %>% 
  select( external_gene_name, feat_fracc_affected, sample, type ) %>% 
  group_by( external_gene_name, type ) %>% 
  summarise( mean_affected = mean( feat_fracc_affected ),
             median_affected = median( feat_fracc_affected ) ) %>% 
  ungroup( )

# join n and aff
summ_all_by_gene.df <- left_join( x = summ_by_gene.df,
                                  y = summ_aff_by_gene.df,
                                  by = c ( "external_gene_name" = "external_gene_name",
                                           "type" = "type" ) )

# save the data for figures
write.table( x = summ_all_by_gene.df, file = "summ_all_by_gene.tsv",
             append = F, quote = F, sep = "\t",
             row.names = F, col.names = T )

# plot 
affected_genes_delcoding.p <- create_scatter_with_marginals_3.f( the_data = summ_all_by_gene.df,
                                                                 filter_type = "protein_coding",
                                                                 the_shape = 25,
                                                                 the_fill = "limegreen",
                                                                 cnv_type = "DEL" )

affected_genes_dellnc.p <- create_scatter_with_marginals_3.f( the_data = summ_all_by_gene.df,
                                                              filter_type = "lncRNA",
                                                              the_shape = 25,
                                                              the_fill = "brown",
                                                              cnv_type = "DEL" )

affected_genes_dupcoding.p <-create_scatter_with_marginals_3.f( the_data = summ_all_by_gene.df,
                                                                filter_type = "protein_coding",
                                                                the_shape = 24,
                                                                the_fill = "limegreen",
                                                                cnv_type = "DUP" )

affected_genes_duplnc.p <- create_scatter_with_marginals_3.f( the_data = summ_all_by_gene.df,
                                                              filter_type = "lncRNA",
                                                              the_shape = 24,
                                                              the_fill = "brown",
                                                              cnv_type = "DUP" )
## save the plot
ggsave( filename = "affected_genes_delcoding.png", plot = affected_genes_delcoding.p, width = 10, height = 10 )
## save the plot
ggsave( filename = "affected_genes_dellnc.png", plot = affected_genes_dellnc.p, width = 10, height = 10 )
## save the plot
ggsave( filename = "affected_genes_dupcoding.png", plot = affected_genes_dupcoding.p, width = 10, height = 10 )
## save the plot
ggsave( filename = "affected_genes_duplnc.png", plot = affected_genes_duplnc.p, width = 10, height = 10 )

# save files required for next module
save( biomart.df, summ_all_by_gene.df, file = "forvenn.RData")
save( ann_affected_genes_FR.df, file = "forfullgenes.RData")
