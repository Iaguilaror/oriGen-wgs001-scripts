/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process venn {

    publishDir "${params.results_dir}/9-3-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
      path( data_channel )
      path( venn_script )

    output:
      path "*", emit: venn_results

    script:
    """
    Rscript --vanilla $venn_script
    """

}

/* name a flow for easy import */
workflow VENN {

 take:
    data_channel
    venn_script

 main:

    venn( data_channel, venn_script )

  emit:
    venn.out[0]

}
