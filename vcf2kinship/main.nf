#!/usr/bin/env nextflow

/*================================================================
The Treviño LAB presents...
  Kinship analysis pipeline
- A tool for finding related samples in VCF files

- This pipeline is meant to reproduce the results in: TO-DO-add url and doi after paper is published

==================================================================
Version: 0.0.1

==================================================================
Authors:
- Bioinformatics Design
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)
 Victor Treviño-Alvarado (vtrevino@tec.mx)

- Bioinformatics Development
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)

=============================
Pipeline Processes In Brief:

Pre-processing:
_000_filter_autosomes

Core-processing:
_001_vcf2plink
_002_king

Pos-processing

Anlysis
_003_plotanalyze

ENDING

================================================================*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE PARAMS DOCUMENTATION AND FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
params.ver = "0.0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
params.pipeline_name = "origen-kinship"

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at JAN 2023
*/
params.nextflow_required_version = '24.04.3'

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.help     = false   //default is false to not trigger help message automatically at every run
params.version  = false   //default is false to not trigger version message automatically at every run

params.input_dir     =	false	//if no inputh path is provided, value is false to provoke the error during the parameter validation block

/* read the module with the param init and check */
include { } from './modules/doc_and_param_check.nf'

// /* load functions for testing env */
include { get_fullParent }  from './modules/useful_functions.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     INPUT PARAMETER VALIDATION BLOCK
  TODO (iaguilar) check the extension of input queries; see getExtension() at https://www.nextflow.io/docs/latest/script.html#check-file-attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = get_fullParent( params.input_dir )

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable (pipeline_name) defined by this Script
  This directories will be automatically created by the pipeline to store files during the run
*/

params.results_dir       =  "${params.output_dir}/${params.pipeline_name}-results/"
params.intermediates_dir =  "${params.output_dir}/${params.pipeline_name}-intermediate/"

/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/* load workflows */
include { KEEPAUTOSOMES }  from  './modules/0-keepautosomes'
include { VCF2PLINK }      from  './modules/1-vcf2plink'
include { KING }           from  './modules/2-king'
include { PLOTKIN }        from  './modules/3-plotanalyze'
include { SELECT2REMOVE }  from  './modules/3b-select-to-remove'
include { REMOVESAMPLES }  from  './modules/4-getsamples-from-vcf'

/* load scripts to send to workdirs */
/* declare scripts channel from modules */
scripts_plotkin = Channel.fromPath( "scripts/king-analysis.R" )
scripts_select2remove = Channel.fromPath( "scripts/king-select.R" )

workflow mainflow {

  main:

// Manually pair the VCF files with their index files
    Channel
    .fromPath( "${params.input_dir}/*.vcf.gz" )
    .map { vcfFile -> [vcfFile, file("${vcfFile}.tbi")] }
	// .view( )
    .set { vcf_channel }

if ( params.run_split_autosomes ){

    king_data = vcf_channel | KEEPAUTOSOMES | VCF2PLINK | KING

} else {

    king_data = vcf_channel | VCF2PLINK | KING

}


    PLOTKIN( king_data, scripts_plotkin )

    samples_to_remove = SELECT2REMOVE( king_data, scripts_select2remove )

if ( params.extract_samples ){
    REMOVESAMPLES( vcf_channel, samples_to_remove )
}

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {

    mainflow( )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/