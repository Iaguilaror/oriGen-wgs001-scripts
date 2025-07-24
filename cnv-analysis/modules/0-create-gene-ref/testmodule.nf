/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { CREATE_GENEREF }    from './main.nf'

 Channel
        .fromPath( "${params.reference}" )
        .set { ref_channel }

/* declare scripts channel for testing */
generef_script = Channel.fromPath( "scripts/0-create-generef.R" )

workflow {
  CREATE_GENEREF( ref_channel, generef_script )
}