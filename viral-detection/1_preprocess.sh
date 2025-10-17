#!/bin/bash

makeblastdb -in test/reference/virus_dumydb.fa -dbtype nucl -out test/reference/virus_dumydb

rm -rf finds/ tmp/

mkdir -p finds
mkdir -p tmp

find test/data -name "*.cram" \
| parallel -j 4 "bash scripts/01_doblast.sh {}" \
&& cat tmp/*.unmap.tmp > finds/unmap_counts.tsv

python generate_summary.py
