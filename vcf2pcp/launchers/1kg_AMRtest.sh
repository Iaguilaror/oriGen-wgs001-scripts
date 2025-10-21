#!/bin/bash

input_dir="test/tutorial-data/"
output_directory="tutorial-results/"
res_config="configfiles/low-res-machine.config"
params_config="configfiles/params-strict.config"

date
nextflow run main.nf \
	--input_dir $input_dir \
	--output_dir $output_directory \
	--sample_list "test/reference/1kg_AMRtest_sample_annotations.tsv" \
	--pca_refpops "test/reference/refpops-for-projection.txt" \
	--projected false \
	--params_config $params_config \
	-c $res_config,$params_config \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html
