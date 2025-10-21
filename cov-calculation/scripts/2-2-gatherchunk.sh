#!/bin/bash

# Get the list of all matching files
files=(*_samtools_mean_dp_all_cov_for_cnv.tsv)

# Extract the index columns
cut -f1-6 "${files[0]}" > base_body.tmp

for f in "${files[@]}"
do

  echo "[DEBUG] pasting $f"
  cut -f7- "$f" > body_inturn.tmp

  paste base_body.tmp body_inturn.tmp > new_base.tmp \
  && mv new_base.tmp base_body.tmp \
  && rm body_inturn.tmp

done

mv base_body.tmp samtools_mean_dp_all_cov_for_cnv.tsv \
&& bgzip samtools_mean_dp_all_cov_for_cnv.tsv
