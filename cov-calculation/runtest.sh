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
	--cram_ref test/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
	--genes_of_interest test/reference/genes_of_interest.txt \
    --contrast_dir test/reference/sample_sets/ \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
