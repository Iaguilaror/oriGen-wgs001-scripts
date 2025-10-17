#!/bin/bash

makeblastdb -in test/reference/virus_dumydb.fa -dbtype nucl -out test/reference/virus_dumydb

rm -rf finds/ tmp/

mkdir -p finds
mkdir -p tmp

find test/data -name "*.cram" | parallel -j 4 "bash scripts/01_doblast.sh {}"

# complete after talking to eugenio

# srun -c 2 python generate_summary.py
