/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process filterlof {

    publishDir "${params.results_dir}/03-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

	tag "$tsv"

    input:
	path(tsv)

    output:
        path "*", emit: filterlof_results

    script:
    """
    awk -F '\t' \
        'NR==1 || (\$34 != "." && \$34 != "")' \
        $tsv > ${tsv.simpleName}_loftee_only.tsv
    """

}

/* name a flow for easy import */
workflow FILTERLOF {

 take:
    tsv_channel

 main:

    tsv_channel | filterlof

  emit:
    filterlof.out[0]

}
