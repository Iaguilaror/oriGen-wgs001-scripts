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
    --seed 100 \
    --admx_threads 1 \
    --amdx_mink 2 \
    --amdx_maxk 3 \
&& echo "[>>>] Module Test Successful"