
/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process gathercovs {

    publishDir "${params.results_dir}/02-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      tuple path( gather_script ), path( mart_channel ), val( chunk_id ), path( counts )

    output:
        path "*", emit: gathercovs_results

    script:
    """
	# Create header with new "sample" column
	echo -e "file\t\$(head -n1 \$(ls *.meandp_byexon.tsv | head -n1) )" > merged_with_sample.tmp

	# Process each file
	for f in *.meandp_byexon.tsv; do
	  sample=\$(basename "\$f" .meandp_byexon.tsv)
	  tail -n +2 "\$f" | awk -v s="\$sample" 'BEGIN{OFS="\t"} {print s, \$0}'
	done >> merged_with_sample.tmp

	bgzip merged_with_sample.tmp

	Rscript --vanilla $gather_script $mart_channel $chunk_id

#	bgzip "$chunk_id"_samtools_mean_dp_all_cov_for_cnv.tsv

	rm *.tmp*

    """

}


process gatherchunks {

	publishDir "${params.results_dir}/02-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

	input:
	  path( chunks )
	  path( chunk_script )

	output:
          path "*", emit: gatherchunks_results

	script:
	"""
	## Rscript --vanilla $chunk_script
	bash $chunk_script
        """
}

/* name a flow for easy import */
workflow GATHERCOVS {

 take:
    tsv_chunks
    gather_script
    mart_channel
    gather_chunks_script

 main:

    all_materials = gather_script
      .combine(mart_channel)
      .combine(tsv_chunks)

    table_chunks = all_materials | gathercovs | toList

    gatherchunks ( table_chunks, gather_chunks_script )

  emit:
    gathercovs.out[0]

}
