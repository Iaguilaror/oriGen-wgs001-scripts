/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { QC }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.tsv" )
        .toList()
        .set { data_channel }

/* declare scripts channel for testing */
qc_script = Channel.fromPath( "scripts/3-QC.R" )

workflow {
  QC( data_channel, qc_script )
}