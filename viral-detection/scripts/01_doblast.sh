#!/bin/bash

the_cram="$1"

mkdir -p finds
mkdir -p tmp
tmpfile=$(mktemp)

python printfasta.py $the_cram "$tmpfile" > tmp/"$(basename $the_cram )".query.fa

blastn \
	-db test/reference/virus_dumydb \
	-query tmp/"$(basename $the_cram )".query.fa -outfmt 6 \
| lz4 -c > finds/"$(basename $the_cram)".tsv.lz4 #\
#> unmap_counts.tsv
