# load packages
pacman::p_load( "vroom", "tidyr", "dplyr",
                "purrr", "stringr" )

# Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "canonical_exons.bed"
# args[2] <- "1"

# put a name to args
biomart <- args[1]
chunk_n <- args[2]

# read the refs
reference.df <- vroom( file = biomart,
                       delim = "\t",
                       col_names = FALSE ) %>% 
  as_tibble( ) %>% 
  rename( chromosome = 1,
          start = 2,
          end = 3,
          gene_name = 4,
          exon_rank = 5,
          transcript_id = 6 )

all_counts.df <- vroom( file = "merged_with_sample.tmp.gz" ) %>%
  unique( ) %>% 
  as_tibble( )

### ----

all_counts_wide.df <- all_counts.df %>%
  pivot_wider(
    id_cols = chrom:id,   # these stay as row identifiers
    names_from = file,              # these become new column names
    values_from = mean_dp           # these fill the cell values
  )

#Debug
# all_counts.df %>%
#   count(gene, exon_rank, file) %>%
#   filter(n > 1)

# fill the DF
final_join.df <- left_join( x = reference.df,
                            y = all_counts_wide.df,
                            by = c( "chromosome" = "chrom",
                                    "start" = "start",
                                    "end" = "end",
                                    "gene_name" = "gene",
                                    "exon_rank" = "exon_rank",
                                    "transcript_id" = "id" ) )

# save the DF
write.table( x = final_join.df,
             file = paste0( chunk_n, "_samtools_mean_dp" , "_all_cov_for_cnv.tsv" ),
             append = F, quote = F, na = "",
             sep = "\t", row.names = F, col.names = T )
