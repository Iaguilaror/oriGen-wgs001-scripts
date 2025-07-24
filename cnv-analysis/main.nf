#!/usr/bin/env nextflow

/*================================================================
The Treviño LAB presents...
  TODO
- A tool for TODO

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
TODO

Core-processing:
TODO

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
params.pipeline_name = "origen-cnv-analysis"

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
params.reference = false

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
include { CREATE_GENEREF }   from  './modules/0-create-gene-ref'
include { CNV_FILTER }       from  './modules/01-filter'
include { GATHER }           from  './modules/02-gather'
include { QC }               from  './modules/03-QC'
include { HEATMAP }          from  './modules/04-heatmap'
include { DOWNLOADMART }     from  './modules/5-biomart-download'
include { BEDTOOLS_CN }      from  './modules/06-bedtools-CN'
include { BEDTOOLS_FR }      from  './modules/07-bedtools-fracc'
//include { GATHER_BED as GATHER_BED_CN } from './modules/08-gather-bedtools'
include { GATHER_BED as GATHER_BED_FR } from './modules/08-gather-bedtools'
include { SUMM_SAMPLE }      from './modules/09-summary-by-sample'


/* load scripts to send to workdirs */
/* declare scripts channel from modules */
generef_script   = Channel.fromPath( "scripts/0-create-generef.R" )
cnvfilter_script = Channel.fromPath( "scripts/1-cnv-filter.R" )
//gather_script    = Channel.fromPath( "scripts/2-gather.R" )
qc_script        = Channel.fromPath( "scripts/3-QC.R" )
heatmap_script   = Channel.fromPath( "scripts/4-heatmap.R" )
biomart_script   = Channel.fromPath( "scripts/0-biomart.R" )
summsample_script= Channel.fromPath( "scripts/9-summsample.R" )
functions_9_script = Channel.fromPath( "scripts/9-functions.R" )

workflow mainflow {

  main:

	Channel
        .fromPath( "${params.reference}" )
        .set { ref_channel }

	ref_channel = CREATE_GENEREF( ref_channel, generef_script )

	Channel
        .fromPath( "${params.input_dir}/*.bcf" )
        .set { bcf_channel }

	data_channel = bcf_channel.combine( cnvfilter_script ).combine( ref_channel )

  base_filter = CNV_FILTER( data_channel )
  .flatten()

	all_tsv = base_filter
  .filter { file -> 
    file.toString( ).endsWith('.filtered_affected_genes_dataframe.tsv') || 
    file.toString( ).endsWith('.filtered_cnv_dataframe.tsv') }
  .toList()

  summary_tsv = GATHER( all_tsv )

  QC( summary_tsv, qc_script )


  cnv_channel = summary_tsv
  .flatten( )
  .filter { file -> file.toString( ).endsWith('allsamples_affected_genes_dataframe.tsv') }

  HEATMAP( cnv_channel, ref_channel, heatmap_script )

  // Do bed analysis
  base_mart = biomart_script | DOWNLOADMART | flatten

  mart_results = base_mart.filter{ it.toString().endsWith('.bed') }
  full_mart = base_mart.filter{ it.toString().endsWith('.tsv') }

  all_bed = base_filter
  .filter { file -> file.toString( ).endsWith('.bed') }

  // add channs
  all_bed
  .combine( mart_results )
  .set{ data_channel }

  all_cn = BEDTOOLS_CN( data_channel ).toList()

  all_fr = BEDTOOLS_FR( data_channel ).toList()

  // Next: gather CN
  //GATHER_BED_CN( all_cn, "allsamples.CN.bed" )
  // Next: gather FR
  gat_FR = GATHER_BED_FR( all_fr, "allsamples.FR.bed" )

  // Next: QC by sample
  // SUMM_SAMPLE( gat_FR, summsample_script, full_mart )
  SUMM_SAMPLE( gat_FR, summsample_script,functions_9_script )
  
  // Next: QC by gene
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
