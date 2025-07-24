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
    --loftee_path /home/bioinfouser/tools/loftee \
    --human_ancestor_fa /home/bioinfouser/tools/loftee_data/human_ancestor_GRCh38_autosomes.fa.gz \
    --conservation_file /home/bioinfouser/tools/loftee_data/loftee.sql \
&& echo "[>>>] Module Test Successful"
