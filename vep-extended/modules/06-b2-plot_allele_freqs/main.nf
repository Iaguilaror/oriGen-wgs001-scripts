/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process gather_af {

    publishDir "${params.results_dir}/06-b-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
	path( materials )

    output:
        path "*", emit: gather_af_results

    script:
    """
    zcat *.allelefrequencies.tsv.gz | awk 'FNR==1 && NR!=1 { next } { print }' > all_chromosomes_AF.tsv
    """

   } 

process af_plot {

    publishDir "${params.results_dir}/06-b-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    tag "$pop"

    input:
	tuple path(all_AF),  path(af_script), val(pop)

    output:
        path "*", emit: af_plot_results

    script:
    """
    Rscript --vanilla $af_script $pop
    """

}

process makepanel {

publishDir "${params.results_dir}/06-c-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    tag "$pop"

    input:
	path( allrds )
    path( panel_script )

    output:
        path "*", emit: makepanel_results

    script:
    """
    Rscript --vanilla $panel_script
    """

}

/* name a flow for easy import */
workflow AF_PLOT {

 take:
    data_channel
    af_script
    panel_script

 main:

    pop_channel = Channel.of(
        'gnomad_AF_afr',
        'gnomad_AF_amr',
        'gnomad_AF_eas',
        'gnomad_AF_nfe',
        'gnomad_AF_sas',
        'mcps_AF_RAW',
        'mcps_AF_IMX'
    )

    all_AF = data_channel | gather_af | combine( af_script ) | combine ( pop_channel )

    allAF_results = af_plot( all_AF ) | flatten

//    allAF_results.view( )

    allrds = allAF_results.filter { it -> it.name.endsWith('.rds') } | toList

    makepanel( allrds , panel_script )

  emit:
    af_plot.out[0]

}
