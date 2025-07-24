/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { VENN }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.RData" )
        .set { data_channel }

/* declare scripts channel for testing */
venn_script = Channel.fromPath( "scripts/9-3-venn.R" )

workflow {
  VENN( data_channel, venn_script )
}