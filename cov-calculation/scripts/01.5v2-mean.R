# load pkgs
pacman::p_load( "vroom", "dplyr" )

# load data
data.df <- vroom( file = list.files( pattern = "\\.sumdp\\.tmp$", full.names = TRUE ),
                  show_col_types = FALSE, col_names = FALSE, na = "." ) %>% 
  as_tibble( ) %>% 
  rename( chrom = 1,
          start = 2,
          end = 3,
          gene = 4,
          exon_rank = 5,
          id = 6,
          read_sum = 7 )

operated.df <- data.df %>% 
  mutate( mean_dp = read_sum / ( end - start - 1 ) ) %>%  # add -1 to emulate result of bedtools map -o mean
  filter( mean_dp > 0 )

# save the data
operated.df %>% 
  # select( gene, exon_rank, mean_dp ) %>% 
  write.table( x = .,
               file = "Rout.tmp",
               append = FALSE, quote = FALSE,
               sep = "\t", row.names = FALSE, col.names = TRUE )
