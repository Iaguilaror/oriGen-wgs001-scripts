/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process countlof {

    publishDir "${params.results_dir}/04-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

	tag "$tsv"

    input:
	    tuple path(tsv), path(script)

    output:
        path "*", emit: countlof_results

    script:
    """
        # create a tmp for lower mem req
        cut -f5,34,38- $tsv \
        | awk '\$2 == "HC" || NR==1' > gt_matrix.tmp

        Rscript --vanilla $script \
            gt_matrix.tmp \
            ${tsv.simpleName}_summary.tsv

        rm gt_matrix.tmp
    """

}

/* name a flow for easy import */
workflow COUNTLOF {

 take:
    tsv_channel

 main:

    tsv_channel | countlof

  emit:
    countlof.out[0]

}