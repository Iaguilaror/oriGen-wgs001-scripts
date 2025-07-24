# load packages
pacman::p_load( "vcfppR", "dplyr", "stringr" )

# load vtrevino code's
source( "https://www.dropbox.com/s/kgszl45a85fqqj4/R2LS.r?dl=1" )

# disable scinot
options( scipen = 6666666 )

# Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "Sample.bcf" # bcf
# args[2] <- "genereference.Rdata" # Rdata ref

# name args
bcf_file <- args[1]
ref_file <- args[2]

# read the ref
load( ref_file )

###################
### DEf functions
getOverlappedGenes <- function (nchr, npos, npos2, naround=10, hgref=hg38) {
  xchr <- paste0("chr",nchr)
  chr <- subset(hgref, chrom == xchr)
  apos <- pmin(npos, npos2)
  bpos <- pmax(npos, npos2)
  chr$dist <- (bpos + apos)/2 - (chr$txEnd + chr$txStart)/2
  chr$OVERLAP <- (!(chr$txEnd < apos | chr$txStart > bpos))
  chrOV <- subset(chr, OVERLAP)
  if (nrow(chrOV)) {
    chrOV$OVSIZE <- pmin(chrOV$txEnd, bpos) - pmax(chrOV$txStart, apos)
    chrOV$OVFRAC <- chrOV$OVSIZE / (chrOV$txEnd - chrOV$txStart)
    o <- order(abs(chrOV$dist))
    return (chrOV[o[1:min(naround,length(o))],c("name","chrom","strand","txStart","txEnd","exonCount","name2","dist","OVERLAP","OVSIZE","OVFRAC")])
  }
  NA
}

annotateOverlappedGenes <- function(chr, pos, pos2, short=FALSE) {
  chr <- sub("^0","",chr)
  a <- character(length(chr))
  for (i in 1:length(a)) {
    x <- getOverlappedGenes(chr[i], pos[i], pos2[i])
    #print(i)
    #print(x)
    if ( typeof( x ) == "list" ) {
      if (short) {
        a[i] <- paste(unique(x$name2[c(1,nrow(x))]),collapse=" ~ ")
      } else {
        a[i] <- paste(unique(x$name2),collapse=":")
      }
    } else {
      a[i] <- ""
    }
  }
  names(a) <- paste0(chr,":",pos)
  a
}

# read the bcffile
###################
### START ANALYSIS

# parametrize for fast allsample_gathering

# BCF - CNV
cnv_bcf <- vcftable( bcf_file, format = "CN" )

# Get svlen from info
cnv_bcf$svlen <- as.numeric( extract( "SVLEN=([^;]+)", cnv_bcf$info, TRUE ) )

# Get svtype from info (DEL, DUP)
cnv_bcf$svtype <- extract( "SVTYPE=([^;]+)", cnv_bcf$info, TRUE )

# Get pytorRD from info (Read Depth Normalizada)
cnv_bcf$rd <- as.numeric( extract( "pytorRD=([^;]+)", cnv_bcf$info, TRUE ) )

# Get pytorP1 from info (e-value (p-val x genomesize / bin size) calcl using t-test between RD statistics in the region and global)
cnv_bcf$p1 <- as.numeric( extract("pytorP1=([^;]+)", cnv_bcf$info,TRUE ) )

# Get pytorP2 from info (e-value (p-val x genomesize / bin size) calcl from the probability of RD values within the region to be in the tails of a gaussian dist of binned RD
cnv_bcf$p2 <- as.numeric( extract("pytorP2=([^;]+)", cnv_bcf$info,TRUE ) )

# Get pytorP3 from info (same as p1 but for middle of CNV)
cnv_bcf$p3 <- as.numeric( extract("pytorP3=([^;]+)", cnv_bcf$info,TRUE ) )

# Get pytorP4 from info (same as p2 but for middle of CNV)
cnv_bcf$p4 <- as.numeric( extract("pytorP4=([^;]+)", cnv_bcf$info,TRUE ) )

# Get pytorQ0 from info (fraq of reads mapped with q0 in call region)
cnv_bcf$q0 <- as.numeric( extract("pytorQ0=([^;]+)", cnv_bcf$info,TRUE ) )

# Get pytorPN from info (fraq of ref genome gaps (Ns) in call region)
cnv_bcf$pN <- as.numeric( extract("pytorPN=([^;]+)", cnv_bcf$info,TRUE ) )

# Get pytorDG from info (distance from closest large (>100bp) gap in reference genome)
cnv_bcf$dG <- as.numeric( extract("pytorDG=([^;]+)", cnv_bcf$info,TRUE ) )

cnv_bcf$end <- as.numeric( extract("END=([^;]+)", cnv_bcf$info, TRUE ) )

# create a dataframe for better handling
cnv_bcf.df <- list.2.data.frame( cnv_bcf )
# gt is the zygosity: 1 HET, 2 HOM

######
##### We should filter here
### From here: https://github.com/abyzovlab/CNVpytor/wiki/2.-Read-Depth-Signal
# "Using viewer mode we can filter calls based on five parameters: CNV size, e-val1, q0, pN and dG:"

# Filter out calls following page info
cnv_bcf_filtered.df <- cnv_bcf.df %>% 
  filter( q0 >= -1 & q0 <= 0.5  ) %>%    ## Filter out calls with more than half not uniquely mapped reads
  filter( p1 >= 0 & p1 <= 0.05  ) %>%    ## filter non-confident calls 
  filter( pN >= 0 & pN <= 0.5  )  %>%    ## filter calls with more than 50% Ns in reference genome  
  # filter( svlen >= 50000 & svlen <= Inf  ) %>%   ## filter calls smaller than 50kbp
  filter( dG >= 100000 & dG <= Inf  )         ##  filter calls close to gaps in reference genome (<100kbp)

# pass the filtered data to the base for syncrony with downstream analysis
cnv_bcf.df <- cnv_bcf_filtered.df

# find genes overlapping the CNV
cnv_bcf.df$ovgenes <- annotateOverlappedGenes( cnv_bcf.df$chr,
                                               cnv_bcf.df$pos,
                                               cnv_bcf.df$pos + cnv_bcf.df$svlen,
                                               short = TRUE )

# find genes overlapping the CNV
cnv_bcf.df$ovgenesall <- annotateOverlappedGenes( chr = cnv_bcf.df$chr,
                                                  pos = cnv_bcf.df$pos,
                                                  pos2 = cnv_bcf.df$pos + cnv_bcf.df$svlen,
                                                  short = FALSE )

# Remove CNVs that affect no genes
cnv_bcf.df2 <- subset( cnv_bcf.df, nchar(ovgenes) > 0 )

# Here save the table of filtered that affect genes
bcf_file %>% 
  str_replace( string = ., pattern = ".bcf$", replacement = ".filtered_cnv_dataframe.tsv" ) %>% 
  write.table( x = select( cnv_bcf.df2, -ref, -qual, -filter, -info),
               file = ., append = F, quote = F,
               sep = "\t", row.names = F, col.names = T  )

# Here save the table of filtered that affect genes, in bed format
bcf_file %>% 
  str_replace( string = ., pattern = ".bcf$", replacement = ".filtered_cnv.bed" ) %>% 
  write.table( x = select( cnv_bcf.df2, chr, pos, end, CN, id, samples ),
               file = ., append = F, quote = F,
               sep = "\t", row.names = F, col.names = F  )

# get a list of all affected genes
sg <- strsplit( cnv_bcf.df2$ovgenesall, ":" )

# create a vector of norm read depth for each affected gene, from the original CNV affecting it
rdg <- rep( cnv_bcf.df2$rd,
            times = unlist( lapply( sg, length ) ) )

# create a datafr with each affected gene and the mean depth of the affecting CNV
dpf <- data.frame( Gene = unlist( sg ),
                   rd = rdg )

# sort by gene name, then by DP
dpf <- dpf[ order( dpf$Gene, dpf$rd ), ]

# sumarize dels & dups by gene
xby <- by( dpf, dpf$Gene, function(x) {
  del <- min( x$rd[ x$rd < 1 ] )
  dup <- max( x$rd[ x$rd > 1 ] )
  data.frame( Gene = x$Gene[ 1 ],
              del = ifelse( test = is.finite( del ),
                            yes = del,
                            no = NA ),
              dup = ifelse( test = is.finite( dup ),
                            yes = dup, no = NA ) )
} )

# prepare for output
odf <- list.2.data.frame( xby, rbind )
# annotate coordinates for affected genes
odf$chrom <- hg38Genes[ odf$Gene, "chrom" ]
odf$start <- hg38Genes[ odf$Gene, "start" ]
odf$end <- hg38Genes[ odf$Gene, "end" ]
odf$sample <- unique( cnv_bcf.df2$samples )

# for output
final.df <- odf %>% 
  select( Gene, del, dup, sample )

# Here save the table of filtered that affect genes
bcf_file %>% 
  str_replace( string = .,
               pattern = ".bcf$",
               replacement = ".filtered_affected_genes_dataframe.tsv" ) %>% 
  write.table( x = final.df,
               file = ., append = F, quote = F,
               sep = "\t", row.names = F, col.names = T  )

# END
