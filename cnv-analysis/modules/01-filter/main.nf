/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process cnv_filter {

    publishDir "${params.results_dir}/1-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      tuple path( bcf ), path( script ), path( ref )

    output:
        path "*", emit: cnv_filter_results

    script:
    """
    Rscript --vanilla $script $bcf $ref
    """

}

/* name a flow for easy import */
workflow CNV_FILTER {

 take:
    data_channel

 main:

    cnv_filter( data_channel )

  emit:
    cnv_filter.out[0]

}
