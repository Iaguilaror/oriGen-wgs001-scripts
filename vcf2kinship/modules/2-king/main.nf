/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process king {

    publishDir "${params.results_dir}/02-king/", mode:"copyNoFollow"

    input:
      tuple path(bed), path(fam), path(bim), path(log)

    output:
        path "king*", emit: king_results

    script:
    """
    king \
      -b $bed \
      --related
    """

}

/* name a flow for easy import */
workflow KING {

 take:
    plink_files

 main:

    plink_files | king 

  emit:
    king.out[0]

}