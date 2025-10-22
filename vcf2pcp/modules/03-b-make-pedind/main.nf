/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process makepedind {

    publishDir "${params.results_dir}/03b-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      path fam_file
      path rscript
      path sample_list

    output:
        path "*", emit: makepedind_results

    script:
    """
    Rscript --vanilla ${rscript} \$(ls *.fam) ${sample_list}
    """

}

/* name a flow for easy import */
workflow MAKEPEDIND {

 take:
    fam_channel
    pedind_script
    sample_list

 main:

    makepedind ( fam_channel, pedind_script, sample_list ) 

  emit:
    makepedind.out[0]

}