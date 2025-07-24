#!/bin/bash

input_dir="real-data/oriGen-p1-1800"
output_directory="real-results/oriGen-p1-1800"
res_config="configfiles/low-res-machine.config"

nextflow run main.nf \
	--input_dir $input_dir \
	--output_dir $output_directory \
	--reference test/reference/hg38-refseq-refgene-ucsc.gz \
	-c $res_config \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html
