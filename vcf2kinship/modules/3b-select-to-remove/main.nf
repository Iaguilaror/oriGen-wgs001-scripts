/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process select2remove {

    publishDir "${params.results_dir}/03b-select2remove/", mode:"copyNoFollow"

    input:
    path MATERIALS

    output:
    path "samples_to_remove.txt", emit: select2remove_results
    path "*.png"

    script:
    """
    Rscript --vanilla king-select.R king.kin
    """

}

/* name a flow for easy import */
workflow SELECT2REMOVE {

 take:

    king_data
    scripts_select2remove

 main:

  king_data
  .combine( scripts_select2remove)
  | select2remove 

  emit:
    select2remove.out[0]

}
