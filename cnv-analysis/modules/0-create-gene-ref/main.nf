/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process create_generef {

    publishDir "${params.results_dir}/0-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      path( ref_channel )
      path( script )

    output:
        path "*", emit: create_generef_results

    script:
    """
    Rscript --vanilla $script $ref_channel
    """

}

/* name a flow for easy import */
workflow CREATE_GENEREF {

 take:
    ref_channel
    generef_script

 main:

    create_generef( ref_channel, generef_script )

  emit:
    create_generef.out[0]

}