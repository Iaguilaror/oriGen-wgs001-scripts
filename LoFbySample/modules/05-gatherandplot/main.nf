/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process gather_plot {

    publishDir "${params.results_dir}/05-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

	tag "$tsv"

    input:
	    path(tsv)
        path(script)

    output:
        path "*", emit: gather_plot_results

    script:
    """
        Rscript --vanilla $script
    """

}

/* name a flow for easy import */
workflow GATHER_PLOT {

 take:
    tsv_channel
    script_gather

 main:

    gather_plot( tsv_channel, script_gather )

  emit:
    gather_plot.out[0]

}