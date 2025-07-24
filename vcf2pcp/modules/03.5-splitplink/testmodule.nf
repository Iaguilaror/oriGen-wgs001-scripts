/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { SPLITPLINK }    from './main.nf'

sample_list =    Channel.fromPath( "${params.sample_list}" )

pcaref_channel = Channel.fromPath( "${params.pca_refpops}" )

    Channel
        .fromPath( "test/data/allwgs.*" )
        .toList( )
        .combine( sample_list )
        .combine( pcaref_channel )
        // .view( )
        .set { data_channel }

/* declare scripts channel for testing */
// NONE

workflow {
  SPLITPLINK( data_channel )
}