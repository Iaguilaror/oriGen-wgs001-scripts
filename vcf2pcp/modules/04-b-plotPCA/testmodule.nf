/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { PLOTPCA }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*" )
        .toList( )
//	.view( )
        .set { data_channel }

/* declare scripts channel for testing */
plotpca_script = Channel.fromPath( "scripts/04b-plotpca.R" )

sample_list = Channel.fromPath( "${params.sample_list}" )

data_channel = data_channel
        .combine( plotpca_script )
        .combine( sample_list )

workflow {
  PLOTPCA( data_channel )
}