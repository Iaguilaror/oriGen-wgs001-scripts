/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { PASTE_ADMX }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*" )
        .toList()
//	.view( )
        .set { pieces_channel }

/* declare scripts channel for testing */
pasteadmx_script = Channel.fromPath( "scripts/05.5-paste_admx.sh" )
pastepopinfo_script = Channel.fromPath( "scripts/05.5-paste_popinfo.R" )

pieces_channel = pieces_channel.combine( pasteadmx_script ).combine( pastepopinfo_script )

workflow {
  PASTE_ADMX( pieces_channel )
}