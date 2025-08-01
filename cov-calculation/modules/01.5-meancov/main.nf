
/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process meancov {

    publishDir "${params.results_dir}/01.5-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      tuple path( sample ), path( bedfile ), path( mean_script )

    output:
        path "*.tsv", emit: meancov_results

    script:
    """
    bedtools map \
      -a $bedfile \
      -b $sample \
      -c 4 \
      -o sum \
      > ${sample.simpleName}.sumdp.tmp

    Rscript --vanilla $mean_script \
    && mv Rout.tmp ${sample.simpleName}.meandp_byexon.tsv

    rm ${sample.simpleName}.sumdp.tmp
    """

}

/* name a flow for easy import */
workflow MEANCOV {

 take:
    data_channel

 main:

    data_channel | meancov 

  emit:
    meancov.out[0]

}
