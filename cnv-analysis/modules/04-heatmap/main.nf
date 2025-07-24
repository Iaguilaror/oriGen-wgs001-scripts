/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process heatmap {

    publishDir "${params.results_dir}/4-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
      path( cnv_channel )
      path( ref_channel )
      path( heatmap_script )

    output:
        path "*", emit: heatmap_results

    script:
    """
    Rscript --vanilla $heatmap_script
    """

}

/* name a flow for easy import */
workflow HEATMAP {

 take:
    cnv_channel
    ref_channel
    heatmap_script

 main:

    heatmap( cnv_channel, ref_channel, heatmap_script )

  emit:
    heatmap.out[0]

}
