/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { PLOTKIN }    from './main.nf'

// Manually pair the bed files with their other files
king_data = Channel.fromPath( "test/data/king.kin" )

/* declare scripts channel for testing */
scripts_plotkin = Channel.fromPath( "scripts/king-analysis.R" )

workflow {
  PLOTKIN( king_data, scripts_plotkin )
}