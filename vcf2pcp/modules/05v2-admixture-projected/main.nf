/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process admixture_projected {

    publishDir "${params.results_dir}/05v2-${task.process.replaceAll(/.*:/, '')}/", mode:"copyNoFollow"

    label 'high_cpu'  // Add a label to the process


    input:
    tuple path(bed), path(bim), path(fam), path(log), path(pfile)

    output:
        path "*", emit: admixture_projected_results

    script:
    """
        K=\$(basename "$pfile" | cut -d. -f2)
        
        cp "$pfile" problem_for_admixture.\$K.P.in

        admixture \
        -j${params.admx_threads} \
        --seed ${params.seed} \
        -P \
        *.bed \
        \$K > problem_admixture.\$K.log
    """

}

/* name a flow for easy import */
workflow ADMIXTURE_PROJECTED {

 take:
    all_channel

 main:

    all_channel
    // .view()
    | admixture_projected

  emit:
    admixture_projected.out[0]

}
