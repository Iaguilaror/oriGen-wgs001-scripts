/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process copyconf {

	publishDir "${params.results_dir}", mode:"copyNoFollow"

 	tag "$conf_channel"  // Add a tag that shows the current instance

    input:
      path conf_channel

    output:
        path "*"

    script:
    """
	echo "# Configs used $conf_channel" > Config_files_and_params.txt
	cat $conf_channel >> Config_files_and_params.txt
    """

}

/* name a flow for easy import */
workflow COPYCONFIG {

 take:
    conf_channel

 main:

    copyconf( conf_channel )

  emit:
    copyconf.out[0]

}
