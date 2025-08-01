/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { GETCOV }    from './main.nf'

/* declare input channel for testing */
// Manually pair the BAM files with their index files
    Channel
        .fromPath( "test/data/*.cram" )
        .map { cramFile -> [cramFile, file("${cramFile}.crai")] }
        .set { cram_channel }

// load reference
Channel
        .fromPath( "test/reference/canonical_exons.bed" )
        .set { bed_channel }

// // load reference
// Channel
//         .fromPath( "${params.cram_ref}" )
//         .set { ref_channel }

// add channs
cram_channel
.combine( bed_channel )
// .combine( ref_channel )
.set{ data_channel }

/* declare scripts channel for testing */
// NONE

workflow {
  GETCOV( data_channel )
}
