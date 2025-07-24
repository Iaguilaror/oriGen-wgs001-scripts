#!/usr/bin/env bash
## This small script runs a module test with the sample data

# remove previous tests
rm -rf .nextflow.log* work

# remove previous results
rm -rf test/results

# create a results dir
mkdir -p test/results

# run nf script
nextflow run testmodule.nf \
	--input_dir test/data/ \
	--dbsnp_ref /bodega/references/dbsnp/156/GRCh38_dbsnp156.vcf.gz \
&& echo "[>>>] Module Test Successful" 
