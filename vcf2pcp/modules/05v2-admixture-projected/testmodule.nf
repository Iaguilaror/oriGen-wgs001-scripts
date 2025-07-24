/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { ADMIXTURE_PROJECTED }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/problem_for_admixture.*" )
        .toList( )
	// .view( )
        .set { data_channel }

    Channel
        .fromPath( "test/reference/*.P" )
        //.toList( )
        .set { p_channel }

    all_channel = data_channel
    .combine( p_channel )
    //.view( )

/* declare scripts channel for testing */
// NONE

workflow {
  ADMIXTURE_PROJECTED( all_channel )
}