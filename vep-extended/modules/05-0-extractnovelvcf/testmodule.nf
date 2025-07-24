/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { EXTRACTNOVEL } from './main.nf'

/* declare input channel for testing */
// Manually pair the VCF files with their index files
    Channel
        .fromPath( "test/data/*.vcf.gz" )
        .set { vcf_channel }

/* declare scripts channel for testing */
// NONE

workflow {
  EXTRACTNOVEL( vcf_channel )
}
