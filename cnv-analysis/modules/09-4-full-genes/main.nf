/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process fullgenes {

    publishDir "${params.results_dir}/9-4-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
      path( data_channel )
      path( fullgenes_script )
      path( functions_9_script )

    output:
      path "*", emit: fullgenes_results

    script:
    """
    Rscript --vanilla $fullgenes_script
    """

}

/* name a flow for easy import */
workflow FULLGENES {

 take:
    data_channel
    fullgenes_script
    functions_9_script

 main:

    fullgenes( data_channel, fullgenes_script, functions_9_script )

  emit:
    fullgenes.out[0]

}
