#!/usr/bin/env nextflow

/*================================================================
The Treviño LAB presents...
  PCA and ADMIXTURE analysis pipeline
- A tool for running PCA and ADMIXTURE in VCF files

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
_001_simplify_removeLD
_002_concatvcf
_003_vcf2plink
_003b_makepedind

Core-processing:
_004_runpca
_004b_plotpca

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
params.pipeline_name = "origen-VCF2PCP-2"

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
params.pca_refpops   =  false

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
include { VCF2PLINK }   from  './modules/03-vcf2plink'
include { MAKEPEDIND }  from  './modules/03-b-make-pedind'

if ( params.projected ) {
  include { SMARTPCA }    from  './modules/04v2-runpca-projected'

  // to "project" admixture we must first split the plink file
  include { SPLITPLINK }  from  './modules/03.5-splitplink/'
  include { ADMIXTURE }   from  './modules/05-admixture/'
  include { ADMIXTURE_PROJECTED }   from  './modules/05v2-admixture-projected/'
  include { PASTE_ADMX }  from  './modules/05.5-paste-projected-admx/'
  include { RENAME_POPINFO }  from  './modules/05.5-b-rename-popinfo/'

} else {
  include { SMARTPCA }    from  './modules/04-runpca'
  include { ADMIXTURE }   from  './modules/05-admixture/'
}

include { PLOTPCA }            from  './modules/04-b-plotPCA'
include { PLOT_ADMIXTURE }     from  './modules/05-b-plot-admixture'
include { GATHER_ADMIXTURE }   from  './modules/05-c-gather-admixture/'

include { COPYCONFIG }         from  './modules/z-copyconfig/'

/* load scripts to send to workdirs */
/* declare scripts channel from modules */
pedind_script       = Channel.fromPath( "scripts/03b-makepedind.R" )
plotpca_script      = Channel.fromPath( "scripts/04b-plotpca.R" )
plotadmx_script     = Channel.fromPath( "scripts/05b-plotadmx.R" )
gat_admx_script     = Channel.fromPath( "scripts/05c-gatheradmx.R" )
pasteadmx_script    = Channel.fromPath( "scripts/05.5-paste_admx.sh" )
pastepopinfo_script = Channel.fromPath( "scripts/05.5-paste_popinfo.R" )

workflow mainflow {

  main:

// Manually pair the VCF files with their index files
    Channel
    .fromPath( "${params.input_dir}/*.vcf.gz" )
    .map { vcfFile -> [vcfFile, file("${vcfFile}.tbi")] }
    .set { vcf_channel }

    /* skip directly to vcf2plink */

    allvcf = vcf_channel

    allplinks = allvcf | VCF2PLINK

    fam_channel = allplinks

    sample_list = Channel.fromPath( "${params.sample_list}" )

    allpedind = MAKEPEDIND( fam_channel, pedind_script, sample_list )


    if ( params.projected ) {

    /* read the poplist */
        Channel
            .fromPath( "${params.pca_refpops}" )
            .set{ pcaref_channel }

        pca_outs = allplinks.combine( allpedind ).combine( pcaref_channel ) | SMARTPCA
    
    /* split to prepare for projected admixture */        

        all_splitplink = allplinks.combine( sample_list ).combine( pcaref_channel ) | SPLITPLINK

        ref_splitplink = all_splitplink.flatten( ).filter { file(it).name.contains('references_for_admixture') }.toList( )
        problem_splitplink = all_splitplink.flatten( ).filter { file(it).name.contains('problem_for_admixture') }.toList( )

    /* run admixture for ref split plink */
        allref_admixresults = ref_splitplink | ADMIXTURE
        p_channel = allref_admixresults.flatten().filter { file(it).name.endsWith('.P') }
        ref_q_channel = allref_admixresults.flatten().filter { file(it).name.endsWith('.Q') }
        
        allprob_admixresults = problem_splitplink
        .combine( p_channel )
        | ADMIXTURE_PROJECTED
        
        prob_q_channel = allprob_admixresults.flatten().filter { file(it).name.endsWith('.Q') }

        pre_popinfo = allpedind
        .flatten( )
        .filter { file(it).name.contains('allwgs.popinfo.txt') }

        pre_fams = all_splitplink
        .flatten( )
        .filter { file(it).name.endsWith('.fam') }

        all_pieces = ref_q_channel
        .mix( prob_q_channel, pre_popinfo, pre_fams, pasteadmx_script, pastepopinfo_script )
        .toList( )

        pre_admixture_outs = all_pieces | PASTE_ADMX

        admixture_outs = pre_admixture_outs
        .flatten()
        .filter { file(it).name.endsWith('.Q') }

        popinfo_channel = pre_admixture_outs
        .flatten()
        .filter { file(it).name.contains('pasted_allwgs.popinfo.txt') }
        | RENAME_POPINFO | flatten


    } else {

        pca_outs = allplinks.combine( allpedind ) | SMARTPCA

        admixture_outs = allplinks | ADMIXTURE

        admixture_outs = admixture_outs
        .flatten()
        .filter { file(it).name.endsWith('.Q') }

        popinfo_channel = allpedind
        .flatten()
        .filter { file(it).name.endsWith('.popinfo.txt') }

    }

    data_channel = pca_outs
        .combine( plotpca_script )
        .combine( sample_list )

    data_channel | PLOTPCA

   all_admix = PLOT_ADMIXTURE( admixture_outs, popinfo_channel, plotadmx_script, sample_list )

   all_admix = all_admix
   .flatten()
   .filter { file(it).name.endsWith('.rds') }
   .toList( )

   GATHER_ADMIXTURE( all_admix, gat_admx_script )

    /* Prepare to copy configs */
  Channel
    .fromPath( "${params.params_config}" )
    .set { conf_channel }

  COPYCONFIG( conf_channel )
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
