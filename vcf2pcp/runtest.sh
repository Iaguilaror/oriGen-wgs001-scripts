#!/bin/bash

input_dir="test/data/"
output_directory="test/results"
res_config="configfiles/low-res-machine.config"
params_config="configfiles/params-strict.config"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run main.nf \
    --input_dir $input_dir \
	--output_dir $output_directory \
	--sample_list "test/reference/OriGen_Phase1_1Kg_and_AiGA.tsv" \
	--pca_refpops "test/reference/refpops-for-projection.txt" \
	--projected true \
	--params_config $params_config \
	-c $res_config,$params_config \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
