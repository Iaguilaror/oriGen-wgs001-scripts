/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process gather {

    publishDir "${params.results_dir}/2-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    input:
      path( data_channel)
      //path( gather_script )

    output:
        path "*", emit: gather_results

    script:
    """
    awk 'FNR==1 && NR!=1 { next } { print }' *_affected_genes_dataframe.tsv > allsamples_affected_genes_dataframe.tsv
    awk 'FNR==1 && NR!=1 { next } { print }' *_cnv_dataframe.tsv > allsamples_cnv_dataframe.tsv

    """

}

/* name a flow for easy import */
workflow GATHER {

 take:
    data_channel
    //gather_script

 main:

    gather( data_channel )

  emit:
    gather.out[0]

}
