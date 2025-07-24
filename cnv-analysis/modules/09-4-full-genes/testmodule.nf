/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { FULLGENES }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.RData" )
        .set { data_channel }

/* declare scripts channel for testing */
fullgenes_script = Channel.fromPath( "scripts/9-4-fullgenes.R" )
functions_9_script = Channel.fromPath( "scripts/9-functions.R" )

workflow {
  FULLGENES( data_channel, fullgenes_script, functions_9_script )
}