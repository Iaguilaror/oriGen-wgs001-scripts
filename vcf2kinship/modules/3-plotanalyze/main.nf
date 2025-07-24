/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process plotking {

    publishDir "${params.results_dir}/03-plotking/", mode:"copyNoFollow"

    input:
    path MATERIALS

    output:
        path "*", emit: plotking_results

    script:
    """
    Rscript --vanilla king-analysis.R king.kin
    """

}

/* name a flow for easy import */
workflow PLOTKIN {

 take:

    king_data
    scripts_plotkin

 main:

  king_data
  .combine( scripts_plotkin)
  | plotking 

  emit:
    plotking.out[0]

}