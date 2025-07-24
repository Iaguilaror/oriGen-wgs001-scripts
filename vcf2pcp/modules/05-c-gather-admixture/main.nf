/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process gather_admixture {

    publishDir "${params.results_dir}/05c-${task.process.replaceAll(/.*:/, '')}/", mode:"copyNoFollow"

    input:
      path MATERIALS
      path RSCRIPT

    output:
        path "*", emit: gather_admixture_results

    script:
    """
    Rscript --vanilla $RSCRIPT
    """

}

/* name a flow for easy import */
workflow GATHER_ADMIXTURE {

 take:
    data_channel
    gat_admx_script

 main:

    gather_admixture( data_channel, gat_admx_script)

  emit:
    gather_admixture.out[0]

}