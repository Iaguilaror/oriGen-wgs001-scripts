/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { MEANCOV }    from './main.nf'

/* declare input channel for testing */
// Manually pair the BAM files with their index files
    Channel
        .fromPath( "test/data/*.readsperposition.bed.gz" )
        .set { data_channel }

// load reference
Channel
        .fromPath( "test/reference/canonical_exons.bed" )
        .set { bed_channel }

/* declare scripts channel for testing */
mean_script = Channel.fromPath( "scripts/01.5v2-mean.R" )

// add channs
data_channel
.combine( bed_channel )
.combine( mean_script )
.set{ data_channel }

/* declare scripts channel for testing */
// NONE

workflow {
  MEANCOV( data_channel )
}
