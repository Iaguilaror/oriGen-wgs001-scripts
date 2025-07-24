# load packages
# pacman::p_load( "vroom", "tidyr", "dplyr" )

# load vtrevino code's
source( "https://www.dropbox.com/s/kgszl45a85fqqj4/R2LS.r?dl=1" )

# Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "hg38-refseq-refgene-ucsc.gz"

# put a name to args
reffile <- args[1]

# load summarized annotation
hg38 <- read.delim( reffile,
                    as.is = TRUE )

hg38Genes <- list.2.data.frame( by( hg38,
                                    hg38$name2,
                                    function( x ) data.frame( Gene = x$name2[ 1 ],
                                                              chrom = x$chrom[1],
                                                              start = min( x$txStart ),
                                                              end = max( x$txEnd ),
                                                              exons = max( x$exonCount ) ) ),
                                rbind )

hg38Genes$size <- hg38Genes$end - hg38Genes$start


# save the hg38 objects
save( hg38, hg38Genes,
      file = "genereference.Rdata")
