/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { FILTERLOF }    from './main.nf'

/* declare input channel for testing */
// Manually pair the VCF files with their index files
    Channel
        .fromPath( "test/data/*.tsv" )
        .set { tsv_channel }

/* declare scripts channel for testing */
// NONE

workflow {
  FILTERLOF( tsv_channel )
}