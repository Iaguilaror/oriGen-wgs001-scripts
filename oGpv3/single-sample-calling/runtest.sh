#!/bin/bash

#WorkflowError in file "/home/iaguilar/Ongoing_projects/oriGen-wgs001-scripts/oGpv3/single-sample-calling/oGpv3.smk", line 8:
#The following environment variables are requested by the workflow but undefined. Please make sure that they are correctly defined before running Snakemake:
#SAMPLES

export SAMPLES="MYRUNID1/MYCODEID1"

snakemake -s oGpv3.smk -c 2
