/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process dbsnp_plot {

    publishDir "${params.results_dir}/06-a-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
	path( materials )

    output:
        path "*", emit: dbsnp_plot_results

    script:
    """
    Rscript --vanilla 06a-dbsnp-analyze.R
    """

}

/* name a flow for easy import */
workflow DBSNP_PLOT {

 take:
    data_channel

 main:

    data_channel | dbsnp_plot

  emit:
    dbsnp_plot.out[0]

}
