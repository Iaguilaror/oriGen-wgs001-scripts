pacman::p_load( "biomaRt", "dplyr", "tidyr" )

# list all marts
listEnsembl( GRCh = 38 )

ensembl <- useEnsembl( biomart = "ensembl" )
all_dts <- listDatasets( ensembl )

# connect to human
ensembl <- useEnsembl( biomart = "ensembl",
                       dataset = "hsapiens_gene_ensembl" )

# list filters
all_filters.df <- listFilters( ensembl )

listFilterOptions( ensembl, "transcript_is_canonical"  )

# useful filters:
# "ensembl_gene_id_version", "ensembl_transcript_id_version",
# "external_gene_name", "biotype", "transcript_biotype",
# "transcript_tsl", "transcript_is_canonical"
# "chromosome_name", 

all_attribts.df <- listAttributes( ensembl )
# useful attrbts:
# feature page
# "ensembl_gene_id_version", "ensembl_transcript_id_version", "chromosome_name",
# "start_position", "end_position", "strand", "transcript_start", "transcript_end",
# "transcript_length", "transcript_tsl", "transcript_is_canonical", "external_gene_name",
#"external_gene_source", "percentage_gene_gc_content", "gene_biotype", "source",
# "external_synonym", "GeneCards ID", "Expression Atlas ID"

# structure
# "ensembl_gene_id_version", "ensembl_transcript_id_version", "ensembl_peptide_id_version",
# "chromosome_name", "start_position", "end_position", "transcript_start",
# "transcript_end", "transcript_length", "strand", "external_gene_name",
# "external_gene_source", "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end",
# "cds_length", "gene_biotype", "exon_chrom_start", "exon_chrom_end", "is_constitutive",
# "rank", "genomic_coding_start", "genomic_coding_end", "ensembl_exon_id, "cds_start", "cds_end"

requested_attr_structure.v <- c(
  "ensembl_transcript_id_version",
  "transcript_version",
  "chromosome_name",
  "start_position",
  "end_position",
  "transcript_start",
  "transcript_end",
  "transcript_length",
  "strand",
  "5_utr_start",                       # structures page
  "5_utr_end",                         # structures page
  "3_utr_start",                       # structures page
  "3_utr_end"                         # structures page
  # "cds_length",
  # "exon_chrom_start",
  # "exon_chrom_end",
  # "is_constitutive",
  # "rank",
  # "genomic_coding_start",
  # "genomic_coding_end",
  # "ensembl_exon_id",
  # "cds_start",
  # "cds_end",
  # "gene_biotype"
)

requested_attr_basic.v <- c(
  "ensembl_gene_id_version",
  "ensembl_transcript_id_version",
  # "ensembl_peptide_id_version",
  # "chromosome_name",
  # "start_position",
  # "end_position",
  # "transcript_start",
  # "transcript_end",
  # "transcript_length",
  "transcript_tsl",
  "transcript_is_canonical",
  "strand",
  # "external_synonym",
  "external_gene_name",
  "external_gene_source",
  "external_transcript_name",
  # "external_transcript_source_name"
  # "cds_length",
  # "exon_chrom_start",
  # "exon_chrom_end",
  # "is_constitutive",
  # "rank",
  # "genomic_coding_start",
  # "genomic_coding_end",
  # "ensembl_exon_id",
  # "cds_start",
  # "cds_end",
  "gene_biotype"
)

filters_list <- list(
  chromosome_name = c(1:22, "X", "Y"),
  transcript_is_canonical = TRUE,
  transcript_tsl = TRUE
)

# basic ts info
transcript.df <- getBM( attributes = requested_attr_basic.v,
                        filters = names( filters_list ),
                        values = filters_list,
                        mart = ensembl )

# count duplicated genes
transcript.df %>% 
  pull( ensembl_transcript_id_version ) %>% 
  duplicated( ) %>% 
  sum( )

# structure ts info
structure.df <- getBM( attributes = requested_attr_structure.v,
                       filters = names( filters_list ),
                       values = filters_list,
                       mart = ensembl )

# get basic location for transcripts
location.df <- structure.df %>% 
  select( ensembl_transcript_id_version:transcript_length ) %>% 
  unique( )

# calculate min max utr # What a mess due to strandness
# Dont do it now
# utrs.df <- structure.df %>% 
#   select( ensembl_transcript_id_version,  )

# join all attri
final_join.df <- left_join( x = location.df,
                           y = transcript.df,
                           by = "ensembl_transcript_id_version" )

# reorder chroms
chrom_order.v <- c(1:22, "X", "Y")

final_join.df <- final_join.df %>% 
  mutate( chromosome_name = factor( chromosome_name, levels = c(1:22, "X", "Y") ) )

# save the full mart for downstream counting
final_join.df %>% 
  write.table( file = "full_biomaRt_used.tsv", append = F, quote = F, sep = "\t", row.names = F, col.names = T )

# save a bedfile for downstream counting
final_join.df %>% 
  select( chromosome_name, start_position, end_position, ensembl_transcript_id_version ) %>% 
  arrange( chromosome_name, start_position, end_position ) %>% 
  write.table( file = "full_biomaRt_used.bed", append = F, quote = F, sep = "\t", row.names = F, col.names = F )