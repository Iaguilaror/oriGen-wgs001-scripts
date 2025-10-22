/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process paste_admx {

    publishDir "${params.results_dir}/05.5-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"


    input:
      path MATERIALS

    output:
        path "*", emit: paste_admx_results

    script:
    
    """
    bash 05.5-paste_admx.sh
    Rscript --vanilla 05.5-paste_popinfo.R
    """

}

/* name a flow for easy import */
workflow PASTE_ADMX {

 take:
    pieces_channel

 main:

    pieces_channel | paste_admx 

  emit:
    paste_admx.out[0]

}