/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { GENEQC }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.bed" )
        .set { bed_channel }

telomers = Channel.fromPath( "${params.telomer}" )
centromers = Channel.fromPath( "${params.centromer}" )

all_references = Channel.fromPath( "test/reference/full_biomaRt_used.tsv" )
.combine( telomers )
.combine( centromers)
.combine( bed_channel ) 

/* declare scripts channel for testing */
geneqc_script = Channel.fromPath( "scripts/9-2-geneqc.R" )
functions_9_script = Channel.fromPath( "scripts/9-functions.R" )

workflow {
  GENEQC( all_references, geneqc_script, functions_9_script )
}