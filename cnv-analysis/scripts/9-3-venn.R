# load packages
pacman::p_load( "vroom", "dplyr", "tidyr", "ggplot2", "stringr", "ggsci",
                "ggExtra", "patchwork", "scales", "pheatmap", "ggvenn" )

# Load objects from prev mod
load( file = "forvenn.RData" )

# get total genes in original mart
all_genes <- biomart.df %>% 
  select( ensembl_gene_id_version, external_gene_name, gene_biotype ) %>% 
  mutate( external_gene_name = ifelse( is.na( external_gene_name ),
                                       yes = ensembl_gene_id_version,
                                       no = external_gene_name ) )

# get DUP genes at least 50%
all_genes_dup.df <- summ_all_by_gene.df %>% 
  filter( type == "DUP" )

# get DEL genes at least 50%
all_genes_del.df <- summ_all_by_gene.df %>% 
  filter( type == "DEL" )

# DO venn for coding

# function to extract genes
extract_gene_names.f <- function( the_data, the_biotype ) {
  
  the_data %>% 
    filter( gene_biotype == the_biotype ) %>% 
    pull( external_gene_name ) %>% 
    unique( )
  
}

all_genes_for_venn_coding.v <- extract_gene_names.f( the_data = all_genes,
                                                     the_biotype = "protein_coding" )


all_genes_for_venn_coding_del.v <- extract_gene_names.f( the_data = all_genes_del.df,
                                                         the_biotype = "protein_coding" )

all_genes_for_venn_coding_dup.v <- extract_gene_names.f( the_data = all_genes_dup.df,
                                                         the_biotype = "protein_coding" )

# propper Venn
# function for venn
do_venn.f <- function( vec1, vec2, vec3, name1, name2, name3 ) {
  
  # Combine the vectors into a named list
  venn_coding.ls <- setNames( list( vec1, vec2, vec3 ), c( name1, name2, name3 ) )
  
  # Create the Venn diagram
  venn_coding.tmp <- ggvenn(
    venn_coding.ls, 
    fill_color = c("tomato", "limegreen", "steelblue"),
    stroke_size = 0.5,
    set_name_size = 5
  ) +
    labs( title = "Gene Set Affected by CNVs",
          caption = "only considering genes affected in at least 50% of its length" ) +
    theme( plot.title =  element_text( hjust = 0.5, size = 15 ) )
  
  return( venn_coding.tmp )
}

venn_coding.p <- do_venn.f( vec1 = all_genes_for_venn_coding_del.v,
                            vec2 = all_genes_for_venn_coding_dup.v,
                            vec3 =  all_genes_for_venn_coding.v,
                            name1 = "OriGEN DEL",
                            name2 = "OriGEN DUP",
                            name3 = "BioMart Protein Coding" )

ggsave( filename = "venn_coding.png", plot = venn_coding.p, width = 10, height = 10, bg = "white" )

# DO venn for lncRNA
all_genes_for_venn_ncRNA.v <- all_genes %>% 
  filter( gene_biotype == "lncRNA" ) %>% 
  pull( external_gene_name ) %>% 
  unique( )

all_genes_for_venn_ncRNA.v <- extract_gene_names.f( the_data = all_genes,
                                                    the_biotype = "lncRNA" )


all_genes_for_venn_ncRNA_del.v <- extract_gene_names.f( the_data = all_genes_del.df,
                                                        the_biotype = "lncRNA" )

all_genes_for_venn_ncRNA_dup.v <- extract_gene_names.f( the_data = all_genes_dup.df,
                                                        the_biotype = "lncRNA" )

venn_ncRNA.p <- do_venn.f( vec1 = all_genes_for_venn_ncRNA_del.v,
                           vec2 = all_genes_for_venn_ncRNA_dup.v,
                           vec3 =  all_genes_for_venn_ncRNA.v,
                           name1 = "OriGEN DEL",
                           name2 = "OriGEN DUP",
                           name3 = "BioMart lncRNA" )

## CHECKPOINT HERE #####
ggsave( filename = "venn_lncRNA.png", plot = venn_ncRNA.p, width = 10, height = 10, bg = "white" )
