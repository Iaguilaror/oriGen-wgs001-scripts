#!/bin/bash

input_dir="real-data/oriGEN-p1-1427"
output_directory="real-results/oriGEN-p1-1427"

res_config="configfiles/low-res-machine.config"

nextflow run main.nf \
	--input_dir $input_dir \
	--output_dir $output_directory \
	--reference test/reference/hg38-refseq-refgene-ucsc.gz \
	-c $res_config \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html
