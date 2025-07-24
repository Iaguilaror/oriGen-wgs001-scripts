# Load packges
pacman::p_load( "purrr", "dplyr", "vroom", "ggplot2", "viridis", "hexbin",
                "ggsci", "scales", "stringr", "tidyr", "ggExtra" )

options( scipen = 666 )

# Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
#args[1] <- "gnomad_AF_amr"

# put a name to args
contrast_pop <- args[1]

alldata.df <- vroom( "all_chromosomes_AF.tsv", na = ".",
                     col_types = c( "c", "n", "c", "c",
                                    rep( "n", 3 ), "c",
                                    rep( "n", 8 )  ) ) %>% 
  as_tibble( ) %>% 
  select( CHROM:oriGen_AN, all_of( contrast_pop ) )

gc( )

# Replace NA with 0 in all numeric columns
alldata.df <- alldata.df %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))

# filter for origen AN
maxan <- max( alldata.df$oriGen_AN )

alldata.df <- alldata.df %>% 
  filter( oriGen_AN >= maxan * 0.9 ) %>% 
  filter( ALT != "*" )

snp_alldata.df <- alldata.df %>% 
  filter( nchar( REF ) == 1 & nchar( ALT ) == 1  ) %>% 
  # select( -CHROM:-ALT, -gnomad_grpmax )
  select( -CHROM:-ALT )

rm( alldata.df )
gc( )

# create funtion to plot
plot_scatter.f <- function( the_df, x_col, y_col, the_caption, the_resolution = 0.005, min_counts = 10 ){
  
  # Define bin edges
  bin_width <- the_resolution
  breaks <- seq( 0, 1, by = bin_width )
  
  # Convert column names to symbols for ggplot
  x_sym <- rlang::sym( x_col )
  y_sym <- rlang::sym( y_col )
  
  # Evaluate symbols to get the column vectors
  x_vec <- rlang::eval_tidy( x_sym, data = the_df )
  y_vec <- rlang::eval_tidy( y_sym, data = the_df )
  
  # Now calculate correlation
  pearson_corr <- cor( x = x_vec,
                       y = y_vec,
                       method = "pearson",
                       use = "complete.obs" )
  
  # Calculate R-squared and format string
  r_squared_label <- paste0( "R\u00B2 = ", round( pearson_corr^2, 3 ) )
  
  # Cut each value into bins
  binned_df <- the_df %>%
    mutate(
      x_bin = cut( !!x_sym, breaks = breaks,
                   include.lowest = TRUE, right = FALSE ),
      y_bin = cut( !!y_sym, breaks = breaks,
                   include.lowest = TRUE, right = FALSE )
    ) %>%
    count( x_bin, y_bin, name = "count" ) %>%
    filter( !is.na( x_bin ) & !is.na( y_bin ) ) %>% 
    filter( count > min_counts )
  
  # Get bin midpoints
  binned_df <- binned_df %>%
    mutate(
      x = as.numeric(sub("\\[(.+),.*", "\\1", x_bin)) + bin_width / 2,
      y = as.numeric(sub("\\[(.+),.*", "\\1", y_bin)) + bin_width / 2
    )
  
  # normalize counts to density = 1
  binned_df <- binned_df %>% 
    mutate( total = sum( count ),
            density = count / total )
  
  density.p <- ggplot( binned_df,
                       mapping = aes(x = x, y = y,
                                     fill = count ) ) +
    geom_tile(  ) +
    geom_abline( intercept = 0,
                 slope = 1,
                 linewidth = 0.5,
                 linetype = "dotted", color = "white" ) +
    scale_fill_gradientn(
      colors =c( "#46deec", "darkblue", "#08306b" ),
      breaks = c(
        # 1e1, #10
        1e2, #100
        # 1e3, #1000
        1e4, #10000
        # 1e5, #100000
        1e6  #1000000
      ),
      labels = c(
        # 1e1, #10
        "100", #100
        # 1e3, #1000
        "10,000", #10000
        # 1e5, #100000
        "1 Million"),  #1000000
        trans = "log10" ) +
    guides(
      fill = guide_colorbar( frame.colour = "black",
                             frame.linewidth = 0.5,   # thinner border,
                             ticks.linewidth = 2,
                             ticks.colour = "white" )
    ) +
    scale_x_continuous( labels = percent ) +
    scale_y_continuous( labels = percent ) +
    labs( x = x_col,
          y = y_col,
          fill = "N Variants",
          title = "2Dhistogram Allele Frequencies in oriGen ~ Other Populations",
          caption = paste0( the_caption, "\nbin_resolution =", the_resolution,
                            "\nmin_counts_to_show =", min_counts ) ) +
    theme_light( base_size = 45 ) +
    theme( panel.grid.minor = element_blank( ),
           panel.background = element_blank( ),
           plot.background = element_blank( ),
           legend.background = element_rect( color = NA, fill = NA ) ) 
  
  p2 <- density.p +
    annotate( "text",
              x = 0.05, y = 0.95,      # Adjust these as needed depending on your data scale
              label = r_squared_label,
              hjust = 0, vjust = 1,
              size = 15 )
  
  ggsave( filename = paste0( the_caption, x_col, "_v_", y_col, "_2dhexbin.svg"),
          plot = p2, width = 15, height = 15 )
  
  p3 <- density.p +
    theme( plot.title = element_blank( ),
           plot.caption = element_blank( ),
           legend.position = "none",
           axis.title = element_blank( ) )
  
  ggsave( filename = paste0( the_caption, x_col, "_v_", y_col, "_2dhexbin_v2.svg"),
          plot = p3, width = 10, height = 10, bg = "transparent", device = svglite::svglite )
  
  # Save a single object to a file
  saveRDS( density.p, paste(contrast_pop, "scatter.rds", sep = "_" ) )
  
  # return the bin to explore
  return( pearson_corr )
  
}

# test func
tmp.df <- snp_alldata.df %>%
  # sample_n( 1e4 ) %>%
  plot_scatter.f( the_df = .,
                  x_col = contrast_pop, y_col =  "oriGen_AF",
                  the_caption = paste( prettyNum( nrow( . ), big.mark = "," ), "SNPs" , sep = "-" ),
                  the_resolution = 0.005, min_counts = 10 )

# # test only density
# tmp.df %>% 
#   filter( density > 0.000001 ) %>% 
# ggplot( data = .,
#         mapping = aes( x = x,
#                        y = y, fill = density ) ) +
#   geom_tile( color = NA )

rm( list = ls( ) )
gc( )
