/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process gather_bed {

    publishDir "${params.results_dir}/X-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
      path( data_channel)
      val( output_file )

    output:
      path "*", emit: gather_bed_results

    script:
    """
    cat *.bed > "$output_file"
    """

}

/* name a flow for easy import */
workflow GATHER_BED {

 take:
    data_channel
    output_file

 main:

    gather_bed( data_channel, output_file )

  emit:
    gather_bed.out[0]

}
