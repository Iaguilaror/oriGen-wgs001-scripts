/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { REMOVESAMPLES }    from './main.nf'

/* declare input channel for testing */
// Manually pair the VCF files with their index files
    Channel
        .fromPath( "test/data/*.vcf.gz" )
        .map { vcfFile -> [vcfFile, file("${vcfFile}.tbi")] }
//	.view( )
        .set { vcf_channel }

/* declare scripts channel for testing */
samples_to_remove = Channel.fromPath( "test/data/samples_to_remove.txt" )

/* declare scripts channel for testing */
// NONE

workflow {
  REMOVESAMPLES( vcf_channel, samples_to_remove )
}