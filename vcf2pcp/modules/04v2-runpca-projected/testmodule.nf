/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { SMARTPCA }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/allwgs.*" )
        .toList( )
//	.view( )
        .set { data_channel }

/* read the poplist */
    Channel
        .fromPath( "${params.pca_refpops}" )
        .set{ pcaref_channel }

/* declare scripts channel for testing */
// NONE

workflow {
  SMARTPCA( data_channel.combine( pcaref_channel ) )
}