/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { GATHER_BED }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.FR.bed" )
        .toList()
        .set { bed_channel }

/* declare scripts channel for testing */
// NONE

workflow {
  GATHER_BED( bed_channel, "allsamples.FR.bed" )
}