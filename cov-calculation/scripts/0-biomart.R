pacman::p_load( "biomaRt", "dplyr", "tidyr" )

# Connect to Ensembl (GRCh38)
ensembl <- useEnsembl( biomart = "ensembl",
                       dataset = "hsapiens_gene_ensembl" )

# Define filters
filters_list <- list(
  chromosome_name = as.character(1:22),      # autosomes only
  transcript_is_canonical = TRUE
)

# Define attributes to retrieve exon structure
requested_attr.v <- c(
  "ensembl_gene_id_version",
  "ensembl_transcript_id_version",
  "external_gene_name",
  "chromosome_name",
  "exon_chrom_start",
  "exon_chrom_end",
  "strand",
  "ensembl_exon_id",
  "rank",
  "is_constitutive",
  "gene_biotype"
)

# Query BioMart
all_exons.df <- getBM(
  attributes = requested_attr.v,
  filters = names(filters_list),
  values = filters_list,
  mart = ensembl
)

# Filter to constitutive exons
all_exons.df <- all_exons.df %>%
  # filter( is_constitutive == 1 ) %>%
  arrange( chromosome_name, exon_chrom_start, exon_chrom_end )

# sort and uniq
final_exons.df <- all_exons.df %>% 
  dplyr::select( chromosome_name, exon_chrom_start, exon_chrom_end, external_gene_name, rank, ensembl_transcript_id_version, gene_biotype ) %>%
  arrange( chromosome_name, exon_chrom_start, exon_chrom_end, external_gene_name, rank ) %>% 
  unique( ) %>% 
#  filter( gene_biotype == "protein_coding" | gene_biotype == "lncRNA" )
  filter( gene_biotype == "protein_coding" ) %>%
  dplyr::select( -gene_biotype )

# Optional: Save as BED
final_exons.df %>%
  write.table( file = "canonical_exons.bed",
               sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
