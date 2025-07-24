# load packages
pacman::p_load( "vroom", "tidyr", "dplyr", "ggplot2", "patchwork", "ggvenn", "scales" )

# load data affected genes
genes.df <- vroom( file = "allsamples_affected_genes_dataframe.tsv" ) %>% 
  as_tibble( )

# get uniq genes
genes_w_del <- genes.df %>% 
  filter( !is.na( del ) ) %>% 
  pull( Gene ) %>% 
  unique( )

samples_w_del <- genes.df %>% 
  filter( !is.na( del ) ) %>% 
  pull( sample ) %>% 
  unique( )

# get uniq genes
genes_w_dup <- genes.df %>% 
  filter( !is.na( dup ) ) %>% 
  pull( Gene ) %>% 
  unique( )

samples_w_dup <- genes.df %>% 
  filter( !is.na( dup ) ) %>% 
  pull( sample ) %>% 
  unique( )

# plot dist of del
del_histo.p <- ggplot( data = genes.df,
                       mapping = aes( x = del ) ) +
  geom_histogram( color = "black", fill = "orange", bins = 100 ) +
  labs( x = "pytorRD (Read Depth Normalized)",
        title = "Depth of DEL",
        subtitle = paste( "From", length( genes_w_del ),
                          "genes, in", length( samples_w_del ), "samples" ) ) +
  theme_classic( base_size = 15 )

dup_histo.p <- ggplot( data = genes.df,
                       mapping = aes( x = dup ) ) +
  geom_histogram( color = "black", fill = "skyblue", bins = 100 ) +
  labs( x = "pytorRD (Read Depth Normalized)",
        title = "Depth of DUP",
        subtitle = paste( "From", length( genes_w_dup ),
                          "genes, in", length( samples_w_dup ), "samples" )) +
  theme_classic( base_size = 15 )

# create a panel
panel_histo.p <- del_histo.p + dup_histo.p

# save the plot
ggsave( filename = "panel_histogram.png",
        plot = panel_histo.p, width = 10, height = 7 )

# Venn diag for affected genes by dup and del
venn_list <- list( DEL = genes_w_del,
                   DUP = genes_w_dup )

venn.p <- ggvenn( venn_list,
                  fill_color = c( "orange", "skyblue" ) ) +
  labs( title = "Number of Unique Affected Genes by CNV type" ) +
  theme( plot.title = element_text( hjust = 0.5, size = 20 ) )

# save the venn
ggsave( filename = "venn_genes.png", bg = "white",
        plot = venn.p, width = 10, height = 7 )

# load data CNV metrics
cnv.df <- vroom( file = "allsamples_cnv_dataframe.tsv" ) %>% 
  as_tibble( )

# plot svlen and mark the recommended cutoff of (50k)
svlen.qc <- ggplot( data = cnv.df,
                    mapping = aes( x = svlen ) ) +
  geom_histogram( bins = 100, color = "black", fill = "white" ) +
  geom_vline( xintercept = 50000,
              lty = "dashed", color = "tomato", alpha = 0.8 ) +
  scale_x_continuous( label = comma ) +
  labs( title = "Range of sizes in CNVs",
        x = "svlen (basepairs)",
        caption = "red line marks cut recommended by CNVpytor
        rec: remove <50k" ) +
  theme_classic( base_size = 12 ) +
  facet_wrap( ~ alt, scales = "free" )

# plot q0 not uniquely mapped reads
q0.qc <- ggplot( data = cnv.df,
                 mapping = aes( x = q0 ) ) +
  geom_histogram( bins = 100, color = "black", fill = "white" ) +
  geom_vline( xintercept = 0.5,
              lty = "dashed", color = "tomato", alpha = 0.8 ) +
  labs( title = "Range of Fracc of not uniquely mapped reads in CNVs",
        x = "q0 ( fracction of not uniquely mapped reads )",
        caption = "red line marks cut recommended by CNVpytor
        rec: remove >0.5" ) +
  theme_classic( base_size = 12 ) +
  facet_wrap( ~ alt, scales = "free" )

# plot pN percentage of Ns in ref genome
pN.qc <-ggplot( data = cnv.df,
                mapping = aes( x = pN ) ) +
  geom_histogram( bins = 100, color = "black", fill = "white" ) +
  geom_vline( xintercept = 0.5,
              lty = "dashed", color = "tomato", alpha = 0.8 ) +
  labs( title = "Range of percentage of Ns in Ref Genome in CNVs",
        x = "pN ( fracction of Ns in ref Genome )",
        caption = "red line marks cut recommended by CNVpytor
        rec: remove >0.5" ) +
  theme_classic( base_size = 12 ) +
  facet_wrap( ~ alt, scales = "free" )

# plot dG distance to gap
dG.qc <- ggplot( data = cnv.df,
                 mapping = aes( x = dG ) ) +
  geom_histogram( bins = 100, color = "black", fill = "white" ) +
  geom_vline( xintercept = 100000,
              lty = "dashed", color = "tomato", alpha = 0.8 ) +
  scale_x_continuous( label = comma ) +
  labs( title = "Range of Distance to Gaps in Genome in CNVs",
        x = "dG ( Distance to gap in ref Genome )",
        caption = "red line marks cut recommended by CNVpytor
        rec: remove <100kbp" ) +
  theme_classic( base_size = 12 ) +
  facet_wrap( ~ alt, scales = "free" )

panel_qc.p <- (svlen.qc + q0.qc) / (pN.qc + dG.qc)

# save the qc
ggsave( filename = "panel_qc.png", bg = "white",
        plot = panel_qc.p, width = 21, height = 10 )

# dg dn qc
dgdn.p <- ggplot( data = cnv.df,
                  mapping = aes( x = pN,
                                 y = log10( dG + 1 ) ) ) +
  geom_point( )

ggsave( filename = "dGdN_qc.png", bg = "white",
        plot = dgdn.p, width = 7, height = 7 )

# create a QC log10
log10dg.p <- ggplot( data = cnv.df,
                     mapping = aes( x = log10( dG + 1 ) ) ) +
  geom_histogram( bins = 100, color = "black", fill = "white" ) +
  geom_vline( xintercept = log10( 100000 + 1 ),
              lty = "dashed", color = "tomato", alpha = 0.8 ) +
  scale_x_continuous( label = comma ) +
  labs( title = "Range of Distance to Gaps in Genome in CNVs",
        x = "log10( dG + 1 ) ( Distance to gap in ref Genome )",
        caption = "red line marks cut recommended by CNVpytor
        rec: remove <100kbp" ) +
  theme_classic( base_size = 12 ) +
  facet_wrap( ~ alt, scales = "free" )

ggsave( filename = "dG_qc.png", bg = "white",
        plot = log10dg.p, width = 10, height = 7 )

# plot CN by type
# Summarize the count of each CN value split by alt
cn_summary_alt <- cnv.df %>%
  group_by(CN, alt) %>%
  summarise(count = n(), .groups = "drop")

# Plot the data
cn_summ.p <- ggplot( data = cn_summary_alt,
                     mapping = aes( x = CN, y = count ) ) +
  geom_bar( stat = "identity",
            position = "dodge",
            color = "black", fill = "white" ) +
  theme_minimal( base_size = 15 ) +
  labs(
    title = "Frequency of Copy Numbers (CN) Split by Type",
    x = "Copy Number (CN)",
    y = "Count",
    fill = "ALT Type"
  ) +
  theme_classic( base_size = 12 ) +
  facet_wrap( ~ alt, scales = "free_y" )

ggsave( filename = "cnsummary_qc.png", bg = "white",
        plot = cn_summ.p, width = 10, height = 7 )
