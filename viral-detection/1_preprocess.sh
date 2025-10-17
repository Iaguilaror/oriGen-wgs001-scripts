#!/bin/bash

makeblastdb -in test/reference/virus_dumydb.fa -dbtype nucl -out test/reference/virus_dumydb

rm -rf finds/ tmp/

mkdir -p finds
mkdir -p tmp

find test/data -name "*.cram" \
| parallel -j 4 "bash scripts/01_doblast.sh {}" \
&& cat tmp/*.unmap.tmp > finds/unmap_counts.tsv

for i in $(ls finds/*.tsv.lz4)
do
    lz4cat $i \
    | perl -lane "print \"$(basename $i .tsv.lz4)\\t\$_\""
done \
| perl scripts/outfmt6_to_bin.pl \
    test/reference/21_Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    test/reference/virus_dumydb.fa.fai \
    test/sampleids.txt

python scripts/generate_summary.py