#!/bin/bash

input_dir="test/data/"
output_directory="test/results"
res_config="configfiles/low-res-machine.config"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run main.nf \
	--input_dir $input_dir \
	--output_dir $output_directory \
	--loftee_path /home/bioinfouser/tools/loftee \
	--human_ancestor_fa /home/bioinfouser/tools/loftee_data/human_ancestor_GRCh38_autosomes.fa.gz \
	--conservation_file /home/bioinfouser/tools/loftee_data/loftee.sql \
	-c $res_config \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
