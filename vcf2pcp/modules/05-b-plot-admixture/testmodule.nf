/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { PLOT_ADMIXTURE }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.Q" )
//      .toList( )
//	.view( )
        .set { data_channel }

/* declare input channel for testing */
    Channel
        .fromPath( "test/aux-data/*.popinfo.txt" )
//      .toList( )
//	.view( )
        .set { aux_channel }

//data_channel = data_channel.combine( aux_channel )

/* declare scripts channel for testing */
plotadmx_script = Channel.fromPath( "scripts/05b-plotadmx.R" )

sample_list = Channel.fromPath( "${params.sample_list}" )

workflow {
  PLOT_ADMIXTURE( data_channel, aux_channel, plotadmx_script, sample_list )
}