/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { DBSNP_PLOT }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.dbsnp_summary.txt" )
        .toList( )
        .set { data_channel }

/* declare scripts channel for testing */
dbsnp_script = Channel.fromPath( "scripts/06a-dbsnp-analyze.R" )

data_channel = data_channel.combine( dbsnp_script )

workflow {
  DBSNP_PLOT( data_channel )
}
