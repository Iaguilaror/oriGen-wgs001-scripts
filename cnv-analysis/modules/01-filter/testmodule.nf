/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { CNV_FILTER }    from './main.nf'

/* declare input channel for testing */
// Manually pair the VCF files with their index files
    Channel
        .fromPath( "${params.input_dir}/*.bcf" )
        // .map { vcfFile -> [vcfFile, file("${vcfFile}.csi")] }
        .set { bcf_channel }

// Read the refence dataframe for genes
    Channel
        .fromPath( "test/reference/genereference.Rdata" )
        .set { ref_channel }

/* declare scripts channel for testing */
cnvfilter_script = Channel.fromPath( "scripts/1-cnv-filter.R" )

data_channel = bcf_channel.combine( cnvfilter_script ).combine( ref_channel )

workflow {
  CNV_FILTER( data_channel )
}