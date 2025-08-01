
/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process plotgenes {

    publishDir "${params.results_dir}/03-${task.process.replaceAll(/.*:/, '')}/", mode:"copyNoFollow"

    input:
      path( data_channel )
      path( gene_script )
      path( goi_channel )
      path( contrast_channel )

    output:
      path "*", emit: plotgenes_results

    script:
    """
    Rscript --vanilla $gene_script $goi_channel
    """

}

/* name a flow for easy import */
workflow PLOTGENES {

 take:
    data_channel
    gene_script
    goi_channel
    contrast_channel

 main:

    plotgenes( data_channel, gene_script, goi_channel, contrast_channel ) 

  emit:
    plotgenes.out[0]

}
