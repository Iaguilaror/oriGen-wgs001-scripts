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
    --gnomad_ref /bodega/references/gnomAD/v4/gnomAD_v4_for_VEP_annotation.vcf.gz \
&& echo "[>>>] Module Test Successful"
