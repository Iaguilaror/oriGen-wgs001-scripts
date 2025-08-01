# load packages
pacman::p_load( "vroom", "tidyr", "dplyr",
                "stringr", "ggplot2" )

# Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "genes_of_interest.txt"

# put a name to args
goi_file <- args[1]

# chose genes of interest
genesof_interest.v <- scan( goi_file, what = character() )

# read contrast file
all_contrast <- list.files( path = ".", pattern = "*.tsv$", full.names = TRUE )

# create a function to read data ans transform
read_transform.f <- function( the_file ){
  
  # read all counts to a single file
  data_of_interest.df <- vroom( file = the_file,
  # data_of_interest.df <- vroom( file = "rpm_all_cov_for_cnv.tsv.gz",
                                show_col_types = FALSE ) %>%
    as_tibble( ) %>%
    filter( gene_name %in% genesof_interest.v )
  
  # select data for plot
  forplot.df <- data_of_interest.df %>% 
    select( gene_name, exon_rank, 7:ncol(.) ) %>%
    mutate( exon_rank = factor( exon_rank, levels = 1:500 ) ) %>% 
    pivot_longer( cols = 3:ncol(.),
                  names_to = "sample",
                  values_to = "mean_dp_samtools"
                  # values_to = "rpm"
                  ) %>% 
    mutate( across( where( is.numeric ), ~ replace_na( ., 0 ) ) ) %>% 
    mutate( gene_name = factor( gene_name, levels = genesof_interest.v ) )
  
  test_samples.df <- forplot.df %>% 
    # filter( sample %in% test_samples.v ) %>% 
    filter( str_detect( sample, str_c( test_samples.v, collapse = "|" ) ) ) %>% 
    unique( ) %>% 
    mutate( tag = "aff_samples" )
  
  ctrl_samples.df <- forplot.df %>% 
    # filter( sample %in% ctrl_samples.v ) %>% 
    filter( str_detect( sample, str_c( ctrl_samples.v, collapse = "|" ) ) ) %>% 
    unique( ) %>% 
    mutate( tag = "ctrl_samples" )
  
  samples_ofinterest.df <- bind_rows( test_samples.df,
                                      ctrl_samples.df )
  
  res_dump.ls <- list( all_samples = forplot.df,
                       test_samples = test_samples.df,
                       ctrl_samples = ctrl_samples.df,
                       interest_samples = samples_ofinterest.df )
  
  return( res_dump.ls )
  
}

# create a function to plot
plot_genes.f <- function( the_list ) {
  
  ggplot( data = the_list$all_samples,
          mapping = aes( x = exon_rank,
                         y = mean_dp_samtools ) ) +
    geom_boxplot( fill = "gray",
                  alpha = 0.5,
                  width = 0.3,
                  outlier.shape = 21, outlier.alpha = 0.3 ) +
    geom_point( data = the_list$interest_samples,
                mapping = aes( x = exon_rank,
                               y = mean_dp_samtools,
                               color = tag ),
                shape = 18, alpha = 0.5,
                position = position_dodge(width = 0.5) ) +
    scale_color_manual( values = c( "tomato", "darkblue" ) ) +
    # labs( caption = "rpm = reads per million exonic reads" ) +
    theme_linedraw() +
    facet_wrap( ~ gene_name, scales = "free_x" ) +
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.y = element_blank() )
  
}

# create a func to replace filenames
oname.f <- function( the_contrast, the_type) {
  
  the_contrast %>% 
    str_remove( pattern = "./" ) %>% 
    str_remove( pattern = ".tsv" ) %>% 
  paste( ., the_type, sep = "_" )
  
}

# this should loop for every contrast asked
for ( contrast_of_interest in all_contrast) {
  
  message( "[DEBUG] running ", contrast_of_interest )
  
  all_samples.df <- vroom( file = contrast_of_interest, show_col_types = FALSE )
  
  test_samples.v <- all_samples.df %>% 
    filter( cnv_condition == "affected" ) %>% 
    pull( sample )
  
  ctrl_samples.v <- all_samples.df %>% 
    filter( cnv_condition == "control" ) %>% 
    pull( sample )
  
  # Plot diff values
  read_transform.f( the_file = "samtools_mean_dp_all_cov_for_cnv.tsv.gz" ) %>%
    plot_genes.f( the_list = . ) %>% 
    ggsave( plot = .,
             filename =  oname.f( the_contrast = contrast_of_interest, the_type = "samtools_mean_dp_boxplot.svg" ),
             width = 7, height = 7 )

  
}
