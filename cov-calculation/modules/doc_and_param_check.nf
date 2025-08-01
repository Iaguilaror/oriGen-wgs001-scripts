/* Define the help message as a function to call when needed *//////////////////////////////
def the_help() {
	log.info"""
  ==========================================
  
  - This pipeline is meant to reproduce the results in: TO-DO-add url and doi after paper is published

  v${params.ver}
  ==========================================
	Usage:
    TO-DO
	""".stripIndent()
}

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	the_help()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "Pipeline v${params.ver}"
	exit 0
}

/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $params.nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $params.nextflow_required_version \n" +
            "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/* Check if inputs provided
    if they were not provided, they keep the 'false' value assigned in the parameter initiation block above and this test fails
*/
if ( !params.input_dir ) {
  log.error " Please provide the following params: --input_dir \n\n" +
  " For more information, execute: nextflow run main.nf --help"
  exit 1
}