/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process qc {

    publishDir "${params.results_dir}/3-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
      path( data_channel)
      path( qc_script )

    output:
        path "*", emit: qc_results

    script:
    """
    Rscript --vanilla $qc_script
    """

}

/* name a flow for easy import */
workflow QC {

 take:
    data_channel
    qc_script

 main:

    qc( data_channel, qc_script )

  emit:
    qc.out[0]

}
