#!/bin/bash

input_dir="real-data/"
output_directory="real-results/1427-ene2025/"

nextflow run main.nf \
	--input_dir $input_dir \
	--output_dir $output_directory \
	--dbsnp_ref /home/bioinfouser/Origen/origen-vepextended/references/GRCh38_dbsnp156_splitmulti.vcf.gz \
	--gnomad_ref /home/bioinfouser/Origen/origen-vepextended/references/gnomAD_v4_for_VEP_annotation.vcf.gz \
	--mcps_ref /home/bioinfouser/Origen/origen-vepextended/references/MCPS_for_VEP_annotation.vcf.gz \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html
