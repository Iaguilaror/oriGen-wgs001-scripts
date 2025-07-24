/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { SUMM_SAMPLE }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.bed" )
        .set { bed_channel }

/* declare scripts channel for testing */
summsample_script = Channel.fromPath( "scripts/9-summsample.R" )
functions_9_script = Channel.fromPath( "scripts/9-functions.R" )

workflow {
  SUMM_SAMPLE( bed_channel, summsample_script, functions_9_script )
}