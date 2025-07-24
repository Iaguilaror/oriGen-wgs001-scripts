#!/bin/bash

input_dir="prepare-env/test/data/"
output_directory="test/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run main.nf \
	--input_dir $input_dir \
	--output_dir $output_directory \
	--dbsnp_ref /bodega/references/dbsnp/156/GRCh38_dbsnp156.vcf.gz \
	--gnomad_ref /bodega/references/gnomAD/v4/gnomAD_v4_for_VEP_annotation.vcf.gz \
	--mcps_ref /bodega/references/mcps/MCPS_for_VEP_annotation.vcf.gz \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
