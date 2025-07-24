/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { AF_PLOT }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.allelefrequencies.tsv.gz" )
        .toList( )
        .set { data_channel }

/* declare scripts channel for testing */
af_script = Channel.fromPath( "scripts/06b-allelefreq-analyze.R" )

panel_script = Channel.fromPath( "scripts/06c-panel.R" )

// data_channel = data_channel.combine( af_script )

workflow {
  AF_PLOT( data_channel, af_script, panel_script )
}
