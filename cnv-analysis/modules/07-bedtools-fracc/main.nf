
/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process bedtools_fr {

    publishDir "${params.results_dir}/07-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      tuple path( samplefile ), path( bedfile )

    output:
        path "*.bed", emit: bedtools_fr_results

    script:
    """
    bedtools intersect \
        -a $bedfile \
        -b $samplefile \
        -wo  \
    | awk '{OFS="\t"; overlap_fraction=\$NF/(\$3-\$2); print \$0, overlap_fraction}' \
    > ${samplefile.simpleName}.FR.bed

    """

}

/* name a flow for easy import */
workflow BEDTOOLS_FR {

 take:
    data_channel

 main:

    data_channel | bedtools_fr

  emit:
    bedtools_fr.out[0]

}