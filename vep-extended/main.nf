#!/usr/bin/env nextflow

/*================================================================
The Treviño LAB presents...
  TBA
- A tool for TBA

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
TBA

Core-processing:
TBA

Pos-processing

Anlysis


ENDING

================================================================*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE PARAMS DOCUMENTATION AND FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
params.pipeline_name = "origen-vcf-freq-annotation"

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     INPUT PARAMETER VALIDATION BLOCK
  TODO (iaguilar) check the extension of input queries; see getExtension() at https://www.nextflow.io/docs/latest/script.html#check-file-attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/* load workflows */
include { CLEANVCF }    from  './modules/00-cleanvcf/'
include { DBSNP }       from  './modules/01-dbSNP/'
include { VEP }         from  './modules/02-vep/'
include { GNOMAD }      from  './modules/03-gnomadv4/'
include { MCPS }        from  './modules/04-mcps/'

include { EXTRACTNOVEL }     from  './modules/05-0-extractnovelvcf/'

include { DBSNP_SUMM }  from  './modules/05-a-dbsnp-summary/'
include { VCF2TSV }     from  './modules/05-b-vcf2tsv/'

include { DBSNP_PLOT }  from  './modules/06-a-plotdbsnp/'
include { AF_PLOT }     from  './modules/06-b-plot_allele_freqs/'

/* load scripts to send to workdirs */
/* declare scripts channel from modules */
dbsnp_script = Channel.fromPath( "scripts/06a-dbsnp-analyze.R" )
af_script = Channel.fromPath( "scripts/06b-allelefreq-analyze.R" )

workflow mainflow {

  main:

// Manually pair the VCF files with their index files
    Channel
    .fromPath( "${params.input_dir}/*.vcf.gz" )
    .map { vcfFile -> [vcfFile, file("${vcfFile}.tbi")] }
    .set { vcf_channel }

    all_clean =  CLEANVCF( vcf_channel )

    all_dbsnp = all_clean | DBSNP

    all_annotated = all_dbsnp | VEP | GNOMAD | MCPS

    all_vcf2tsv = all_annotated | VCF2TSV

    all_dbsnp_summ = all_dbsnp | DBSNP_SUMM

    all_dbsnp_summ.toList( ).combine( dbsnp_script ) | DBSNP_PLOT

    tsv_channel = all_vcf2tsv.toList( )

    AF_PLOT( tsv_channel, af_script )

    // extract novel VCF
    EXTRACTNOVEL( all_annotated )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {

    mainflow( )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
