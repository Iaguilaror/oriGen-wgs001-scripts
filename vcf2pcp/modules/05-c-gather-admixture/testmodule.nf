/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { GATHER_ADMIXTURE }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.rds" )
        .toList( )
//	.view( )
        .set { data_channel }

/* declare scripts channel for testing */
gat_admx_script = Channel.fromPath( "scripts/05c-gatheradmx.R" )

workflow {
  GATHER_ADMIXTURE( data_channel, gat_admx_script )
}