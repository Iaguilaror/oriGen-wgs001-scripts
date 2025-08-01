/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { GATHERCOVS }    from './main.nf'

/* declare input channel for testing */
    Channel
        .fromPath( "test/data/*.meandp_byexon.tsv" )
        .toSortedList( )
        .set { tsv_channel }

collated_tsv = tsv_channel
   .flatten( )
   .collate( params.table_chunk_size )

group_ids = collated_tsv
    .toList( )
    .map { it.size() }
    .map { N -> (1..N).toList() }     // creates list [1, 2, ..., N]
    .flatten()                        // emits each number as a separate item

/*
collated_tsv
//    .toList()
    .view()

group_ids
    .view()
*/

tsv_chunks = group_ids
   .merge( collated_tsv ) { a, b -> tuple(a, b) }
//    .view( )

/* get the mart reference to extract biotypes */
mart_channel =Channel
        .fromPath( "test/reference/canonical_exons.bed" )

/* declare scripts channel for testing */
gather_script = Channel
	.fromPath( "scripts/2-gather.R" )

gather_chunks_script = Channel
//        .fromPath( "scripts/2-2-gatherchunk.R" )
	.fromPath( "scripts/2-2-gatherchunk.sh" )

workflow {
  GATHERCOVS( tsv_chunks, gather_script, mart_channel, gather_chunks_script )
}

