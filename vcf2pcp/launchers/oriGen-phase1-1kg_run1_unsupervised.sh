#!/bin/bash

input_dir="real-data/origen-p1-run1_1427_1kg_NOmissingto0_MAF0.05_noLD0.8/"
output_directory="real-results/origen-p1-run1_1427_1kg_NOmissingto0_MAF0.05_noLD0.8-projected-500random_origen_on_training_uptoK7_notprojectedPCA/"
res_config="configfiles/high-res-machine.config"
params_config="configfiles/params-strict.config"

date
nextflow run main.nf \
	--input_dir $input_dir \
	--output_dir $output_directory \
	--sample_list "test/reference/origen-p1-1427-1kg.tsv" \
	--pca_refpops "test/reference/refpops-for-projection.txt" \
	--projected false \
	--params_config $params_config \
	-c $res_config,$params_config \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html
