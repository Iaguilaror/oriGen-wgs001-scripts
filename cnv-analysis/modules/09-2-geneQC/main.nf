/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process geneqc {

    publishDir "${params.results_dir}/9-2-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
      path( all_materials )
      path( geneqc_script )
      path( functions_9_script )

    output:
      path "*", emit: geneqc_results

    script:
    """
    Rscript --vanilla $geneqc_script
    """

}

/* name a flow for easy import */
workflow GENEQC {

 take:
    all_materials
    geneqc_script
    functions_9_script

 main:

    geneqc( all_materials, geneqc_script, functions_9_script )

  emit:
    geneqc.out[0]

}
