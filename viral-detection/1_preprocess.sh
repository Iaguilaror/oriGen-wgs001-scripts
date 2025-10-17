#!/bin/bash

makeblastdb -in test/reference/virus_dumydb.fa -dbtype nucl -out test/reference/virus_dumydb


# srun -c 2 python generate_summary.py
