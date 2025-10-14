#!/usr/bin/env nextflow

/*================================================================
The Treviño LAB presents...
  Coverage comparisson analysis pipeline
- A tool for calculating coverage on CRAMS from different samples and compare

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
_000_biomart

Core-processing:
_001_getcov
_002_gather

Pos-processing

Anlysis


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
params.pipeline_name = "origen-cov-calculation"

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at SEP 2024
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
include { DOWNLOADMART }   from  './modules/0-biomart-download'

include { GETCOV }         from  './modules/01-getcov'
include { MEANCOV }         from  './modules/01.5-meancov'

include { GATHERCOVS }     from  './modules/02-gather_covs'

/* load scripts to send to workdirs */
/* declare scripts channel from modules */
biomart_script  = Channel.fromPath( "scripts/0-biomart.R" )
mean_script     = Channel.fromPath( "scripts/01.5v2-mean.R" )
gather_script   = Channel.fromPath( "scripts/2-gather.R" )

gather_chunks_script = Channel.fromPath( "scripts/2-2-gatherchunk.sh" )

gene_script     = Channel.fromPath( "scripts/3-geneplot.R" )

workflow mainflow {

  main:

// Manually pair the CRAM files with their index files
    Channel
    .fromPath( "${params.input_dir}/*.cram" )
    .map { cramFile -> [cramFile, file("${cramFile}.crai")] }
    .set { cram_channel }

    mart_channel = DOWNLOADMART( biomart_script )

	// add channs
	cram_channel
	.combine( mart_channel )
	.set{ data_channel }

  cov_channel = GETCOV( data_channel )

  // add channs
  cov_channel
  .combine( mart_channel )
  .combine( mean_script )
  .set{ data_channel_2 }

  all_mean = MEANCOV( data_channel_2 )  | toList

  collated_tsv = all_mean
   .flatten( )
   .collate( params.table_chunk_size )

  group_ids = collated_tsv
    .toList( )
    .map { it.size() }
    .map { N -> (1..N).toList() }     // creates list [1, 2, ..., N]
    .flatten()

tsv_chunks = group_ids
   .merge( collated_tsv ) { a, b -> tuple(a, b) }


  tsv_channel = GATHERCOVS( tsv_chunks, gather_script, mart_channel, gather_chunks_script )

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
