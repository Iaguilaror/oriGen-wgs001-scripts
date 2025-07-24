/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { GATHER_PLOT }    from './main.nf'

/* declare input channel for testing */
// Manually pair the VCF files with their index files
    Channel
        .fromPath( "test/data/*.tsv" )
        .toList( )
        .set { tsv_channel }

/* declare scripts channel for testing */
script_gather = Channel
    .fromPath( "scripts/05-gather.R" )

workflow {
  GATHER_PLOT( tsv_channel, script_gather )
}
