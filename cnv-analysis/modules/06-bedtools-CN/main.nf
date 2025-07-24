
/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process bedtools_cn {

    publishDir "${params.results_dir}/06-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      tuple path( samplefile ), path( bedfile )

    output:
        path "*.bed", emit: bedtools_cn_results

    script:
    """
    bedtools intersect \
        -a $bedfile \
        -b $samplefile \
        -wa -wb  \
    > ${samplefile.simpleName}.CN.bed

    """

}

/* name a flow for easy import */
workflow BEDTOOLS_CN {

 take:
    data_channel

 main:

    data_channel | bedtools_cn 

  emit:
    bedtools_cn.out[0]

}