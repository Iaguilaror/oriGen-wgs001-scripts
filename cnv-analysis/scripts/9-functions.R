# To gather all functions used to explore results of CNV

###
# from 9-summsample.R ####
# Function to create a column plot for total affected BP per sample
create_bars_ordered_bysample.f <- function(the_data, col_to_analyze,
                                           the_type = "DEL",
                                           the_fill = "gray", 
                                           the_alpha = 1,
                                           the_title = "Total Affected BP per Sample", 
                                           the_x_title = "Sample",
                                           the_y_title = "Total Affected BP") {
  col_to_analyze <- enquo( col_to_analyze )
  
  # Filter the dataframe and prepare the data for plotting
  plot_data <- the_data %>%
    filter(type == the_type) %>%              # Filter based on the_type
    arrange(!!col_to_analyze) %>%             # Sort by the column to analyze
    mutate(sample = factor(sample, levels = sample)) # Convert sample to a factor for sorting
  
  # Create the column plot
  the_bars <- ggplot(plot_data, aes(x = sample, y = !!col_to_analyze)) +
    geom_col(fill = the_fill, alpha = the_alpha) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x-axis labels
    labs(
      title = the_title,
      x = the_x_title,
      y = the_y_title
    ) +
    scale_y_continuous(
      labels = comma,
      sec.axis = sec_axis(~ ., name = the_y_title, labels = comma) # Duplicate the y-axis with the same labels
    ) +
    scale_x_discrete( expand = c( 0.02, 0.02 ) ) +
    theme_linedraw(base_size = 15) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  
  return( the_bars )
}

# Function to create scatterplot with marginal boxplots
create_scatter_with_marginals.f <- function( data, x_col, y_col,
                                             the_fill = "gray", the_shape = 21,  the_alpha = 1,
                                             the_size = 3, base_size = 15,
                                             the_title, the_subtitle, x_title, y_title, the_caption = "",
                                             the_hline = NA ) {
  
  # Create the scatterplot
  scatter_plot <- ggplot(data = data, mapping = aes_string(x = x_col, y = y_col)) +
    geom_point(shape = the_shape, fill = the_fill, size = the_size, alpha = the_alpha ) +
    theme_minimal( ) +
    labs(
      title = the_title,
      subtitle = the_subtitle,
      x = x_title,
      y = y_title,
      caption = the_caption
    ) +
    theme_linedraw(base_size = base_size)
  
  if( !is.na( the_hline ) ){
    scatter_plot <- scatter_plot +
      geom_hline( yintercept = the_hline, alpha = 0.5,
                  lty = "dashed", color = "tomato", size = 1 )
  }
  
  # Add marginal boxplots
  scatter_with_marginals <- ggMarginal(
    scatter_plot,       # The base ggplot object
    type = "boxplot",   # Type of marginal plot
    margins = "both",   # Marginals for both x and y axes
    size = 10,          # Size of marginal plots
    alpha = the_alpha,
    fill = the_fill,  # Fill color for the boxplots
    alpha = 0.5
  )
  
  return(scatter_with_marginals)
}

# summary for numbers
summ_ngenes.f <- function( the_data, the_biotype, the_cnvtype ) {
  the_data %>% 
    filter( gene_biotype == the_biotype ) %>%
    filter( type == the_cnvtype ) %>% 
    summary( )
}

#
# Function to create scatterplot with marginal boxplots
create_scatter_with_marginals_3.f <- function( the_data,
                                               filter_type = "protein_coding",
                                               the_fill = "gray",
                                               the_shape,
                                               cnv_type ) {
  
  filtered_data <- the_data %>%
    filter( gene_biotype == filter_type) %>%
    filter( type == cnv_type)
  
  # Create the base scatterplot
  scatter_plot <- ggplot( data = filtered_data,
                          mapping = aes( x = n, y = mean_affected ) ) +
    geom_vline( xintercept = 1427 * 0.01, color = "tomato", alpha = 0.5, size = 2, lty = "dashed" ) +
    geom_vline( xintercept = 1427 * 0.05, color = "orange", alpha = 0.5, size = 2, lty = "dashed" ) +
    geom_point( alpha = 0.2, shape = the_shape, fill = "black" ) +
    scale_x_log10( breaks = c(10, 100, 1000, 10000),  # Custom breaks for x-axis
                   labels = scales::label_number( accuracy = 1, big.mark = "" ) ) + # Labels as plain numbers
    scale_y_continuous( labels = scales::percent ) + # Percent format for y-axis
    labs(
      title = "Dist of Genes by CNV alteration",
      subtitle = paste( cnv_type, "-", filter_type ),
      x = "Number of samples (log10 scale)",
      y = "Mean % of affected gene length",
      caption = "Only considering genes affected in at least 50% of its length"
    ) +
    theme_linedraw( base_size = 15 )
  
  # Add marginal boxplots
  scatter_with_marginals <- ggMarginal(
    scatter_plot,
    type = "boxplot",
    margins = "both",
    size = 10,
    fill = the_fill,
    alpha = 0.5
  )
  
  return( scatter_with_marginals )
  
}

# From 9-4 Fullgenes

# to transfor fill value
unlog10.f <- function( the_log ){ 10^the_log %>% floor }

# function to plot 
do_rangeplot.f <- function( the_data, the_biotype ) {
  
  gene_data.tmp <- the_data %>% 
    filter( gene_biotype == the_biotype ) %>% 
    arrange( min_CN, delta ) %>% 
    mutate( external_gene_name = factor( external_gene_name, levels = external_gene_name  ) )
  
  ggplot( data = gene_data.tmp,
          mapping = aes( x = external_gene_name ) ) + 
    geom_segment( mapping = aes( xend = external_gene_name, y = min_CN, yend = max_CN ),  alpha = 0.1 ) +
    geom_point( mapping = aes( y = mean_CN, fill = log10( samples ) ), shape = 21, alpha = 0.5 ) +
    scale_y_continuous( limits = c( 0, 26 ),
                        breaks = seq( 0, 26, 2 ) ) +
    scale_x_discrete( expand = c( 0.05, 0.05 ) ) +
    scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 2,
                          limits = c( 1, 3.2 ),
                          breaks = c( 1, 2, 3 ),
                          labels = unlog10.f,
                          guide = guide_colorbar( 
                            frame.colour = "black", 
                            frame.linewidth = 0.5,
                            ticks.colour = "black",  # Add black ticks
                            ticks.linewidth = 0.5      # Set tick line width
                          ) ) +
    labs( title = paste( the_biotype, "- Range of Copy Number in Whole Genome" ),
          subtitle = "Dot shows Mean CN, gray lines show Min to Max CN",
          caption = "Only in 100% affected genes\nOnly genes in >1%",
          x = "Gene",
          y = "CN",
          fill = "N Samples" ) +
    theme_classic( base_size = 15 ) +
    theme( axis.text.x = element_blank( ) )
  
}

# function to plot a bar venn diagram
barvenn.f <- function( the_data, the_subtitle ){
  # set 1: DUPs
  fullgenes_duplicated.v <- the_data %>% 
    filter( min_CN >= 3 | max_CN >= 3  ) %>% 
    pull( external_gene_name ) 
  
  # set 2: DELs
  fullgenes_deleted.v <- the_data %>% 
    filter( min_CN <= 1 | max_CN <= 1  ) %>% 
    pull( external_gene_name ) 
  
  # set 0: ALL
  fullgenes_all.v <- c( fullgenes_duplicated.v, fullgenes_deleted.v ) %>% 
    unique( )
  
  # set 3: common
  # Find common elements
  fullgenes_common.v <- intersect(fullgenes_duplicated.v, fullgenes_deleted.v)
  
  # create a bar venn diagram
  for_bar_venn.df <- tibble( "set" = c( "Duplicated Genes", "Deleted Genes", "Dup OR Del Genes" ),
                             "N_genes" = c( length(fullgenes_duplicated.v), length( fullgenes_deleted.v ), length( fullgenes_common.v ) ) ) %>% 
    mutate( min_y = case_when( set == "Duplicated Genes" ~ 0,
                               set == "Deleted Genes" ~ length( fullgenes_all.v) - N_genes,
                               set == "Dup OR Del Genes" ~ length(fullgenes_duplicated.v) - N_genes ) ) %>% 
    mutate( max_y = case_when( set == "Duplicated Genes" ~ N_genes,
                               set == "Deleted Genes" ~ length( fullgenes_all.v),
                               set == "Dup OR Del Genes" ~ min_y + N_genes ) ) %>% 
    mutate( min_x = 2,
            max_x = 3 ) %>% 
    unique( ) %>% 
    mutate( mid_y = case_when( set == "Duplicated Genes" ~ max( min_y ) / 2 ,
                               set == "Deleted Genes" ~ min( max_y ) + ((max( max_y ) - min( max_y )) / 2) ,
                               set == "Dup OR Del Genes" ~ max( min_y ) + (N_genes / 2) ) ) %>% 
    mutate( text = case_when( set == "Duplicated Genes" ~ N_genes - length( fullgenes_common.v ),
                              set == "Deleted Genes" ~ N_genes - length( fullgenes_common.v ),
                              set == "Dup OR Del Genes" ~ N_genes ) ) %>% 
    mutate( the_label = paste( str_remove( string = set, pattern = " Genes" ),
                               text, sep = "\n") )
  
  # Create a rectangle plot
  my_bar_colors = c( "Duplicated Genes" = "blue", "Deleted Genes" = "red", "Dup OR Del Genes" = "black" )
  
  the_plot <- ggplot( data = for_bar_venn.df,
                        mapping = aes(
                          xmin = min_x,
                          xmax = max_x,
                          ymin = min_y,
                          ymax = max_y,
                          fill = set ) ) +
    geom_rect(  size = 0.5, alpha = 0.5) + # Black border and semi-transparent fill
    geom_text( mapping = aes( x = 3.5,
                              y = mid_y,
                              label = the_label ),
               size = 8 ) +
    scale_y_continuous( name = "Number of Genes" ) + # Y-axis label
    scale_x_continuous( limits = c( 0.5, 4.5 ) ) + # X-axis label
    scale_fill_manual( values = my_bar_colors ) +
    labs( title = "Genes affected by CNVs",
          x = "Gene Set",
          caption = "For genes 100% DELeted or DUPlicated",
          subtitle = the_subtitle ) +
    coord_flip( ) +
    theme_minimal( base_size = 15 ) +
    theme( legend.position = "none",
           axis.line.x = element_line( ),
           axis.ticks.x = element_line( ),
           axis.title.y = element_text( angle = 0, vjust = 0.5 ),
           panel.grid = element_blank( ),
           axis.text.y = element_blank( ) ) # Optional: position the legend at the top
  
  the_df <- bind_rows( data.frame( genes =  fullgenes_common.v,
                                    condition = "DUP OR DEL" ),
                        data.frame( genes = setdiff( fullgenes_duplicated.v,
                                                     fullgenes_deleted.v ),
                                    condition = "DUP exclusive" ),
                        data.frame( genes = setdiff( fullgenes_deleted.v,
                                                     fullgenes_duplicated.v ),
                                    condition = "DEL exclusive" )  )
  
  results.list <- list( the_plot = the_plot,
                        the_df = the_df )
}

# heatmaps
make_pheatmap.f <- function( the_data, the_biotype ){
  
  # Make the hatmap for CNs
  full_cnv_matrix.df <- the_data %>% 
    filter( gene_biotype == the_biotype ) %>% 
    select( sample, external_gene_name, CN ) %>% 
    unique( ) %>% 
    pivot_wider( id_cols = external_gene_name,
                 names_from = "sample",
                 values_from = "CN",
                 values_fill = 2 )
  
  # Convert the tibble to a matrix for pheatmap (excluding the external_gene_name column)
  cnv_matrix <- as.matrix(full_cnv_matrix.df[ , -1])
  
  # Set row names to the external_gene_name column
  rownames(cnv_matrix) <- full_cnv_matrix.df$external_gene_name
  
  # Generate the heatmap
  pheatmap(cnv_matrix, 
           color = colorRampPalette(c("white", "orange", "red"))(50), 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           show_rownames = FALSE, 
           show_colnames = FALSE,
           main = paste( "CopyNumber for fully (100%) DUP or DEL genes", "-", the_biotype ) )
  
  # Add custom axis labels
  # grid::grid.text("Genes", x = 0.05, y = 0.5,
  #                 rot = 90, gp = grid::gpar(fontsize = 14))
  # grid::grid.text("Samples", x = 0.5, y = 0.9,
  #                 rot = 0, gp = grid::gpar(fontsize = 14))
  
}
