/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process plot_admixture {

    publishDir "${params.results_dir}/05b-${task.process.replaceAll(/.*:/, '')}/", mode:"copyNoFollow"

    input:
      path MATERIALS

    output:
        path "*", emit: plot_admixture_results

    script:
    """
    Rscript --vanilla $MATERIALS
    """

}

/* name a flow for easy import */
workflow PLOT_ADMIXTURE {

 take:
    data_channel
    aux_channel
    plotadmx_script
    sample_list

 main:

    plotadmx_script.combine( data_channel ).combine( aux_channel ).combine( sample_list ) | plot_admixture

  emit:
    plot_admixture.out[0]

}