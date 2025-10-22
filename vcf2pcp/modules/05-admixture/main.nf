/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process admixture {

    publishDir "${params.results_dir}/05-${task.process.replaceAll(/.*:/, '')}/", mode:"copyNoFollow"

    tag "K=$K"  // Add a tag that shows the current K value

    label 'high_cpu'  // Add a label to the process

//    maxForks 1 // Limit to N concurrent instances

    input:
//      path MATERIALS
//     val K
    tuple path(bed), path(bim), path(fam), path(log), val(K)

    output:
        path "*", emit: admixture_results

    script:
    """
      admixture \
        -j${params.admx_threads} \
        --seed ${params.seed} \
        *.bed \
        $K > admixture.'$K'.log
    """

}

/* name a flow for easy import */
workflow ADMIXTURE {

 take:
    data_channel

 main:

//    all_k = Channel.from( 4..params.amdx_maxk.toInteger() )
    all_k = Channel.from( params.amdx_mink.toInteger()..params.amdx_maxk.toInteger() )

    data_channel
    .combine(all_k)
    // .view()
    | admixture

  emit:
    admixture.out[0]

}
