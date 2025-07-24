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
    --telomer test/reference/telomere_positions.txt \
    --centromer test/reference/centromeres.tsv \
&& echo "[>>>] Module Test Successful" 