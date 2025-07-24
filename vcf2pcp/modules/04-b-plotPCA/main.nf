/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process plotpca {

    if ( params.projected ) {
      out_prefix = "04b-projected"
    } else {
      out_prefix = "04b-not-projected"
    }

    publishDir "${params.results_dir}/${out_prefix}-${task.process.replaceAll(/.*:/, '')}/", mode:"copyNoFollow"

    input:
      path MATERIALS

    output:
        path "*", emit: plotpca_results

    script:
    """
    Rscript --vanilla 04b-plotpca.R \
    "allwgs.evec.gz" \
    "allwgs.eval" \
    "allwgs.tracy_widom_statistics" \
    "smartpca.stdout" \
    "parallel_plot.svg" \
    "\$(basename ${params.sample_list})"
    """

}

/* name a flow for easy import */
workflow PLOTPCA {

 take:
    data_channel

 main:

    data_channel | plotpca 

  emit:
    plotpca.out[0]

}