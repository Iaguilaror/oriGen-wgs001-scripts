#!/bin/bash

input_dir="real-data/ten_samples"
output_directory="real-results/ten_samples"
res_config="configfiles/low-res-machine.config"
# params_config="configfiles/params-strict.config"

nextflow run main.nf \
    --input_dir $input_dir \
	--output_dir $output_directory \
	--cram_ref test/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
	--genes_of_interest test/reference/genes_of_interest.txt \
	--contrast_dir test/reference/sample_sets/ \
	-c $res_config \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html
