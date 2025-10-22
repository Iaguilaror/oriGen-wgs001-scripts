#!/bin/bash

input_dir="test/tutorial-data/"
output_directory="test/tutorial-results/"
res_config="configfiles/low-res-machine.config"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run main.nf \
    --input_dir $input_dir \
	--output_dir $output_directory \
	--sample_list "test/tutorial-reference/1kg_AMRtest_sample_annotations.tsv" \
	--pca_refpops "test/tutorial-reference/refpops-for-projection.txt" \
	--projected true \
	--amdx_mink 4 \
	--amdx_maxk 6 \
	--admx_threads 2 \
	--seed 100 \
	-c $res_config \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
