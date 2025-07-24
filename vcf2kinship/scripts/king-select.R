# load packages
pacman::p_load( "vroom", "dplyr",
                "ggplot2", "ggsci",
                "patchwork", "tidyverse",
                "scales" )

# Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "test/data/king.kin"

# put a name to args
ifile <- args[1]

# read data
king_data <- vroom( file = ifile )

# select data
king_data_select <- king_data %>% 
  select( ID1, ID2, InfType ) %>% 
  arrange( InfType )

###
# focus on close related idns
close_realted.df <- king_data_select %>% 
  filter( InfType == "Dup/MZ" | InfType == "PO" | InfType == "FS" | InfType == "2nd" )

close_realted.df2 <- close_realted.df %>% 
  mutate( id_1 = ID2,
          id_2 = ID1 ) %>% 
  mutate( ID1 = id_1,
          ID2 = id_2 ) %>% 
  select( -id_1, -id_2 )

close_fullmatrix <- bind_rows( close_realted.df,
                               close_realted.df2 )

### Find samples to remove
# create a function to count samples in pairs
sample.counts.f <- function( the_dataframe ) {
  
  # Count the number of appearances of each sample
  counts <- the_dataframe %>%
    pivot_longer(cols = c(ID1, ID2), names_to = "ID_type", values_to = "Sample") %>%
    count(Sample) %>%
    arrange(desc(n))
  
  return( counts )
  
}

# Get the original sample list full
original_samples <- sample.counts.f( the_dataframe = king_data_select ) %>% 
  pull( Sample )

# count
sample_counts <- sample.counts.f( the_dataframe = close_fullmatrix )
  
remaining_pairs <- close_fullmatrix

# init a vector for removed samples
removed_samples <- vector( "character" )

# Function to iteratively remove samples
while ( nrow(remaining_pairs) > 0) {
  
  # Identify the sample with the highest count
  sample_to_remove <- sample_counts$Sample[1]
  
  # send debug mess
  message( "[DEBUG] removing ", sample_to_remove )
  
  # Remove all pairs involving this sample
  remaining_pairs <- remaining_pairs %>%
    filter( !( ID1 == sample_to_remove | ID2 == sample_to_remove ) )
  
  # Update sample counts
  sample_counts <- sample.counts.f( the_dataframe = remaining_pairs )
  
  # register the sample to remove
  removed_samples <- c( removed_samples, sample_to_remove )
  
  if ( nrow( sample_counts ) == 0 ) break
  
}

# compare before and after removing samples ----
king_data_original_counts <- king_data_select %>% 
  group_by( InfType ) %>% 
  summarise( n = n( ) )

king_data_unrelated_counts <- king_data_select %>% 
  filter( ! (ID1 %in% removed_samples | ID2%in% removed_samples) ) %>% 
  group_by( InfType ) %>% 
  summarise( n = n( ) )

# all samples
keep_samples <- original_samples[ ! original_samples %in% removed_samples ]

# save the sample names
write.table( x = sort( keep_samples ), file = "samples_to_keep.txt",
             append = F, quote = F,
             sep = "\t" , row.names = F, col.names = F )

write.table( x = sort( removed_samples ), file = "samples_to_remove.txt",
             append = F, quote = F,
             sep = "\t" , row.names = F, col.names = F )

# create a simple plot for a bar of eliminated samples
plot.df <- data.frame( type = c( "unrelated", "related"),
                       samples = c( length( keep_samples ), length( removed_samples ) ),
                       axis = c( "project", "project" ) ) %>% 
  mutate( percentage = samples / sum( samples ) ) %>% 
  mutate( label = paste( type, 
                         prettyNum( samples, big.mark = "," ),
                         percent( percentage ), sep = " | " ) ) %>% 
  arrange( type )

# myscale
myscale.f <- function( the_number ){
  prettyNum( the_number, big.mark = "," )
}

# plot data
barra_related.p <- ggplot( data = plot.df,
        mapping = aes( x = axis,
                       y = samples,
                       fill = type ) ) +
  geom_col( color = "black", width = 0.2 ) +
  scale_fill_manual( limits = plot.df$type,
                     labels = plot.df$label,
                     name = "Kinships",
                     values = c( "related" = "tomato",
                                 "unrelated" = "skyblue" ) ) +
  scale_y_continuous( labels = myscale.f ) +
  labs( title = "Related samples in the dataset",
        caption = "Reported by king v2" ) +
  theme_classic( ) +
  theme( axis.title.x = element_blank( ) )

# save plot
ggsave( filename = "barra_related.png", plot = barra_related.p,
        width = 5, height = 7 )
