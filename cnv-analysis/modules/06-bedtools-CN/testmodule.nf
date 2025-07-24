/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { BEDTOOLS_CN }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.bed" )
        .set { sample_channel }

// load reference
Channel
        .fromPath( "test/reference/full_biomaRt_used.bed" )
        .set { bed_channel }

// add channs
sample_channel
.combine( bed_channel )
.set{ data_channel }

/* declare scripts channel for testing */
// NONE

workflow {
  BEDTOOLS_CN( data_channel )
}