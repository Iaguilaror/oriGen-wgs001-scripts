# oriGen Phase 1 – Reproducible Analysis Code

- This pipelines are meant to reproduce results from: TO-DO-add url and doi after paper is published

Reproducible code and workflows for the oriGen Phase 1 analysis of 1,427 whole genomes.  

Welcome! 👋  
This repository contains code used to analyze the data in the oriGen Phase 1 study, which involved whole-genome sequencing of 1,427 individuals. Our goal is to provide transparent and reproducible access to the workflows behind our results.

## 📄 About the Study

The oriGen Phase 1 study is a population-scale genomic project focused on Mexican individuals. This release includes:

- Variant Effect Annotation
- CNV Analysis
- Loss of Function Summary
- Kinship Analysis
- ADMIXTURE and Principal Component Analysis

## How to use this repository

Each workflow subdirectory contains the Nextflow code used to run analysis on a part of the paper. No sensitive nor personal data is shared in the repositories.

Every subdirectory includes a ./scripts/ dir with particular R, python or bash scripts for descriptive and or statistical analyses.

Every subdirectory includes a main.nf, and a modules/ subdir with the nextflow framework to run.

## 📁 Workflow Subdirectory Structure

```bash
├── vep-extended/    # Variant Effect Annotation
├── cnv-analysis/    # CNV Analysis
├── LoFbySample/     # Loss of Function Summary
├── vcf2kinship/     # Kinship Analysis
├── vcf2pcp/         # ADMIXTURE and Principal Component Analysis
└── README.md        # This file
```
