/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process summ_sample {

    publishDir "${params.results_dir}/9-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
      path( bed_channel )
      path( summsample_script )
      path( functions_9_script )

    output:
        path "*", emit: summ_sample_results

    script:
    """
    Rscript --vanilla $summsample_script
    """

}

/* name a flow for easy import */
workflow SUMM_SAMPLE {

 take:
    bed_channel
    summsample_script
    functions_9_script

 main:

    summ_sample( bed_channel, summsample_script, functions_9_script )

  emit:
    summ_sample.out[0]

}
