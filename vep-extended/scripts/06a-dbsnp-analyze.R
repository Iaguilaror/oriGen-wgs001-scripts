# Load packges
pacman::p_load( "purrr", "dplyr", "vroom", "ggplot2",
                "ggalluvial", "ggsci", "scales", "stringr" )

# find all files
allsumaries <- list.files( path = ".",
                           pattern = "*.dbsnp_summary.txt",
                           full.names = TRUE )

# load all files
alldata.df <- allsumaries %>% 
  map_df( ~ vroom( ., col_types = c("c", "c", "n", "c" ) ) ) 

# plot a bar for known novel
novel_summ <- alldata.df %>% 
  group_by( dbsnp, type ) %>% 
  summarise( n = sum( n_variants ) ) %>% 
  ungroup( ) %>% 
  filter( dbsnp != "total" ) %>% 
  group_by( dbsnp ) %>% 
  mutate( dbsnp_total = sum( n )  ) %>% 
  ungroup(  ) %>% 
  mutate( percentage = dbsnp_total / ( sum( dbsnp_total ) / 2 ) ) %>% 
  mutate( dbsnp = paste( dbsnp,
                         percent( percentage, accuracy = 0.1 ),
                         sep = "\n" ) ) %>% 
  mutate( type = str_replace( string = type, pattern = "snp", replacement = "SNP" ) ) %>% 
  mutate( type = str_replace( string = type, pattern = "indel", replacement = "INDEL" ) )

forcaption <- novel_summ %>% 
  select( type, dbsnp, n ) %>% 
  mutate( dbsnp = str_replace( string = dbsnp,
                              pattern = "\n", replacement = " " ) )
  

# Create the alluvial plot
sankey.p <- ggplot( data = novel_summ,
        mapping = aes( axis1 = dbsnp,
                       axis2 = type,
                       y = n ) ) +
  geom_alluvium( mapping = aes( fill = dbsnp ),
                 width = 1/3,
                 knot.pos = 1/4,
                 knot.prop = TRUE,
                 curve_type = NULL,
                 curve_range = NULL,
                 segments = NULL,
                 na.rm = FALSE,
                 show.legend = NA,
                 inherit.aes = TRUE ) +
  geom_stratum( mapping = aes( fill = dbsnp ), # <--- map fill here too
                width = 0.2, alpha = 0.5 ) +
  geom_text( stat = "stratum",
             mapping = aes( label = after_stat( stratum ) ) ) +
  scale_fill_d3( ) +
  labs( title = "OriGen Whole Genomes",
        caption = "compared to dbsnp 156        " ) +
  theme_void( ) +
  theme( plot.title = element_text( hjust = 0.5 ),
         plot.subtitle = element_text( hjust = 0.5 ),
         legend.position = "none" )

# save plot and dataframe
ggsave( filename = "sankey_novel.png",
        plot = sankey.p, bg = "white",
        width = 5, height = 5 )

# save df
write.table( x = forcaption, file = "novel_variants.tsv",
             append = F, quote = F, sep = "\t",
             row.names = F, col.names = T )
