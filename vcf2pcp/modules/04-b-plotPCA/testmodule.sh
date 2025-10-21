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
    --sample_list "test/reference/1kg_AMRtest_sample_annotations.tsv" \
    --projected true \
&& echo "[>>>] Module Test Successful"
