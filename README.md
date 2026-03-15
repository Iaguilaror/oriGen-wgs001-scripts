# oriGen – Analysis Code

- The code in this repository are meant to reproduce results from: TO-DO-add url and doi after paper is published

Workflows for the oriGen Analysis of 1,427 whole genomes.  

Welcome! 👋  
This repository contains code used to analyze the data in the oriGen study, which involved whole-genome sequencing of 1,427 individuals. Our goal is to provide transparent and reproducible access to the workflows behind our results.

## 📄 Contents of this repository  

The oriGen study is a population-scale genomic project focused on Mexican individuals. This release includes the following subdirectories:  

- **oGpv3**: NGS preprocessing, QC, alignment and Variant Calling
- **vep-extended**: Variant Effect Annotation
- **LoFbySample**: Loss of Function variation summary
- **vcf2kinship**: Kinship Analysis with KING2
- **vcf2pcp**: ADMIXTURE and Principal Component Analysis from joint population VCF
- **cnv-analysis**: Processing of the CNV BCFs from CNVpytor
- **cov-calculation**: NGS coverage calculation by exon for CN estimation
- **viral-detection**: Detection of viral hits in sequencing data

## How to use this repository

Each workflow subdirectory contains:

* The Source code used to run analysis on a part of the paper
* A small (simulated or real) dataset to demo the code
* All software dependencies and operating systems  
* Installation instructions
* Output descriptions
* Expected runtimes

### * **No sensitive nor personal data is shared in the repositories** *

Every subdirectory includes a ./scripts/ dir with particular R, python or bash scripts for descriptive and or statistical analyses.

Every subdirectory includes a main.nf, and a modules/ subdir with the nextflow framework to run; or a snakemake equivalent.

## 📁 Workflow Subdirectory Structure

```bash
├── cnv-analysis/    # CNV Analysis from CNVpytor output
├── cov-calculation  # Coverage calculation for CNV estimation
├── LoFbySample/     # Loss of Function Summary
├── oGpv3            # pipeline for pre-processing, QC, alignment and short variant calling
├── vcf2kinship/     # Kinship Analysis
├── vcf2pcp/         # ADMIXTURE and Principal Component Analysis
├── vep-extended/    # Variant Effect Annotation
├── viral-detection  # Analysis of viral hits in NGS data
├── files_for_projection_ADMXandPCA # Variant list, ADMIXTURE P files, and PCA weights to project your samples
└── README.md        # This file
```

### Contact
If you have questions, requests, or bugs to report, please open an issue in this github page.  

### Dev Team
Israel Aguilar-Ordonez <iaguilaror@gmail.com>   
Victor Trevino Alvarado <vtrevino@tec.mx>   
Eugenio Guzman Cerezo <eugenio.guzman@tec.mx>   

This code was developed as part of Israel Aguilar-Ordoñez’s postdoctoral research at Tecnológico de Monterrey during the 2024–2025.

### Cite us

If you find the code in this repository useful, please include the following citation in your work:

TO-DO after publication.  
