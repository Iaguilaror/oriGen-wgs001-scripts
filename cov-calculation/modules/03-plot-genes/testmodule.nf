/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { PLOTGENES }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.tsv.gz" )
        .toList( )
        .set { data_channel }

Channel
        .fromPath( "${params.genes_of_interest}" )
        .set { goi_channel }

Channel
        .fromPath( "${params.contrast_dir}/*.tsv" )
        .toList( )
        .set { contrast_channel }

/* declare scripts channel for testing */
gene_script = Channel.fromPath( "scripts/3-geneplot.R" )

workflow {
  PLOTGENES( data_channel, gene_script, goi_channel, contrast_channel )
}