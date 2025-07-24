
/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process downloadmart {

    publishDir "${params.results_dir}/5-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      path( biomart_script )

    output:
        path "*", emit: downloadmart_results

    script:
    """
    Rscript --vanilla $biomart_script

    """

}

/* name a flow for easy import */
workflow DOWNLOADMART {

 take:
    biomart_script

 main:

    biomart_script | downloadmart 

  emit:
    downloadmart.out[0]

}