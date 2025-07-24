# load packages
pacman::p_load( "vroom", "tidyr", "dplyr" )

# Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "allwgs.converted2plink.fam"
# args[2] <- "OriGen_Phase1_1Kg_and_AiGA.tsv"


# put a name to args
ifile <- args[1]
samples <- args[2]

# read data
fam.df <- vroom( file = ifile, col_names = FALSE )
sample.df <- vroom( file = samples, comment = "#" )

# join data
merged.df <- fam.df %>% 
  left_join( x = .,
             y = sample.df,
             by = c( "X2" = "sample" ) ) %>% 
  mutate( `Superpopulation code` = ifelse( test = is.na( `Superpopulation code` ),
                                           yes = "UNKNOWN",
                                           no = `Superpopulation code` ) )

tosave.df <- merged.df %>% 
  separate( col = X1,
            into = c( "project", "sample" ),
            sep = "-", remove = FALSE ) %>% 
  mutate( X1 = paste( project, X2, sep = "-" ) ) %>% 
  select( X1:X5, `Superpopulation code`, -project, -sample )

# Display the updated dataframe
head(fam.df)


# save the data
write.table( x = tosave.df, file = "allwgs.pedind", 
             append = FALSE, quote = FALSE, 
             sep = " ", row.names = FALSE, col.names = FALSE )

# save the data
tosave.df %>%
  select( 1,2,6 ) %>%
  write.table( x = ., file = "allwgs.popinfo.txt", 
               append = FALSE, quote = FALSE, 
               sep = "\t", row.names = FALSE, col.names = FALSE )
