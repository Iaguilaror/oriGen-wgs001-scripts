/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process rename_popinfo {

    publishDir "${params.results_dir}/05.5-b-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"


    input:
      path MATERIALS

    output:
        path "*", emit: rename_popinfo_results

    script:
    
    """
    cp pasted_allwgs.popinfo.txt allwgs.popinfo.txt
    """

}

/* name a flow for easy import */
workflow RENAME_POPINFO {

 take:
    data_channel

 main:

    data_channel | rename_popinfo 

  emit:
    rename_popinfo.out[0]

}