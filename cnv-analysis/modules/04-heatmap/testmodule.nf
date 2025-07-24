/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { HEATMAP }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/allsamples_affected_genes_dataframe.tsv" )
        .toList()
        .set { cnv_channel }

// Read the refence dataframe for genes
    Channel
        .fromPath( "test/reference/hg38-refseq-refgene-ucsc.gz" )
        .set { ref_channel }

/* declare scripts channel for testing */
heatmap_script = Channel.fromPath( "scripts/4-heatmap.R" )

workflow {
  HEATMAP( cnv_channel, ref_channel, heatmap_script )
}