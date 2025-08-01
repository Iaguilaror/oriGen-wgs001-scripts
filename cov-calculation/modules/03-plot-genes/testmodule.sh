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
    --genes_of_interest test/reference/genes_of_interest.txt \
    --contrast_dir test/reference/sample_sets/ \
&& echo "[>>>] Module Test Successful" 