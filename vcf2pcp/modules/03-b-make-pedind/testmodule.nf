/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { MAKEPEDIND }    from './main.nf'

/* declare input channel for testing */
// Manually pair the VCF files with their index files
    Channel
        .fromPath( "test/data/*.fam" )
//	.view( )
        .set { fam_channel }

/* declare scripts channel for testing */
pedind_script = Channel.fromPath( "scripts/03b-makepedind.R" )

sample_list = Channel.fromPath( "${params.sample_list}" )

workflow {
  MAKEPEDIND( fam_channel, pedind_script, sample_list )
}