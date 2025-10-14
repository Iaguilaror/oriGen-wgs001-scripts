# cnv-analysis
Repo to analyze Copy Number Variants pre-calculated with CNVpytor

===  

- A tool for analyzing CNVpytor sample level results  

- This pipeline is meant to reproduce the results in: TO-DO-add url and doi after paper is published

'cnv-analysis' is a pipeline tool that takes copy number variant data in BCF format to generate the following outputs:  

````
For MULTIPLE BCF and .csi pair of files

1) allsamples_affected_genes_dataframe.tsv  # Summary of genes affected by sample. Long format
2) allsamples_cnv_dataframe.tsv             # All CNVs detected, by sample, annotated with genes affecteds
3) DEL_heatmap.png                          # Visual clustering for affected genes, all samples. For Deletions
4) DUP_heatmap.png                          # Visual clustering for affected genes, all samples. For Duplications
5) bp_delbysample.png                       # Distribution of deleted base pairs
6) bp_dupbysample.png                       # Distribution of duplicated base pairs
7) bp_scatter.png                           # Scatterplot of DUP and DEL affected base pairs
8) n_CNV_by_sample.tsv                      # Summary of the number of CNVs in each sample
9) n_cnv_scatterside.png                    # Scatterplot of number of DUP and DEL in all samples
10) n_delbysample.png                       # Barplots for number of deletions by sample
11) n_dupbysample.png                       # Barplots for number of duplications by sample

````
---

### Features
  **-v 0.0.1**

* Supports BCF and .csi files
* Results include tables with All CNVs by sample, annotated with affected genes
* Results include a plots for QC
* Scalability and reproducibility via a Nextflow-based framework   

---

## Requirements
#### This pipeline was successfully teste on the following OS*:  
* [Ubuntu 22.04.4 LTS](https://releases.ubuntu.com/focal/)

#### Incompatible OS*:
* UNKNOWN  

\* origen-LoFbySample may run in other UNIX based OS and versions, but testing is required.  

#### Command line Software required:
| Requirement | Version  | Required Commands * |
|:---------:|:--------:|:-------------------:|
| [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 24.04.3 | nextflow |
| [bedtools](https://anaconda.org/bioconda/bedtools) | 2.31.1 | bedtools |
| [R](https://www.r-project.org/) | 4.4.1 (2024-06-14) | Rscript |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).  

#### R packages required:

```
vroom version: 1.6.5
tidyr version: 1.3.1
dplyr version: 1.1.4
stringr version: 1.5.1
ggplot2 version: 3.5.1
patchwork version: 1.3.0
scales version: 1.4.0
ggsci version: 3.2.0
ggExtra version: 0.10.1
biomaRt version: 2.62.1
vcfppR version: 0.7.6
ggvenn version: 0.1.10
pheatmap version: 1.0.12
tibble version: 3.2.1
```

---

### Installation
Download pipeline from Github repository:  
```
git clone https://github.com/Iaguilaror/oriGen-wgs001-scripts.git

cd oriGen-wgs001-scripts/cnv-analysis
```

---

## Replicate our analysis (Testing the pipeline):

* Estimated test time:  **3 minute(s)**  
* on a 16 core, 64G RAM machine  

1. To test pipeline execution using test data, run:  
```
bash runtest.sh
```

2. Your console should print the Nextflow log for the run, once every process has been submitted, the following message will appear:  
```
======
 Basic pipeline TEST SUCCESSFUL
======
```

3. Pipeline results for test data should be in the following directory:  
```
./test/results/
```
---


### Pipeline Inputs

* A directory containing a `.bcf file and .vcf.gz.tbi index` with genotypes of multiple samples.

Example contents  
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=2023-11-16
##reference=GRCh38
##source=CNVpytor
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
...
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  samplecnv
1       2       CNVpytor_del0   N       <DEL>   .       PASS    END=10000;IMPRECISE;SVLEN=10000;SVTYPE=DEL;pytorRD=0;pytorP1=0;pytorP2=0;pytorP3=3.98431e-11;pytorP4=0;pytorQ0=-1;pytorPN=1;pytorDG=0;pytorCL=rd_mean_shift     GT:CN   1/1:0
1       50001   CNVpytor_del1   N       <DEL>   .       PASS    END=90000;IMPRECISE;SVLEN=40000;SVTYPE=DEL;pytorRD=0.661052;pytorP1=7.96863e-12;pytorP2=5.67113e-09;pytorP3=2.51641e-11;pytorP4=4.35399e-08;pytorQ0=0.408956;pytorPN=0;pytorDG=40100;pytorCL=rd_mean_shift      GT:CN   0/1:1
...
```  

---

### Pipeline Results

```
Inside the directory test/results/origen-vepextended/2-gather you can find the following:  

1) allsamples_affected_genes_dataframe.tsv  # Summary of genes affected by sample. Long format
2) allsamples_cnv_dataframe.tsv             # All CNVs detected, by sample, annotated with genes affecteds

Inside the directory test/results/origen-vepextended/4-heatmap you can find the following:  

3) DEL_heatmap.png                          # Visual clustering for affected genes, all samples. For Deletions
4) DUP_heatmap.png                          # Visual clustering for affected genes, all samples. For Duplications

Inside the directory test/results/origen-vepextended/9-summ_sample you can find the following:  

5) bp_delbysample.png                       # Distribution of deleted base pairs
6) bp_dupbysample.png                       # Distribution of duplicated base pairs
7) bp_scatter.png                           # Scatterplot of DUP and DEL affected base pairs
8) n_CNV_by_sample.tsv                      # Summary of the number of CNVs in each sample
9) n_cnv_scatterside.png                    # Scatterplot of number of DUP and DEL in all samples
10) n_delbysample.png                       # Barplots for number of deletions by sample
11) n_dupbysample.png                       # Barplots for number of duplications by sample

```


---

### module directory structure

````
.
├── main.nf         # the Nextflow main script
├── modules/        # sub-dirs for development of the Nextflow modules
├── README.md       # This readme
├── runtest.sh      # bash script to launch the pipeline test locally
├── scripts/        # directory with all the scripts used by the pipeline
└── test
    └── data       # sample data to run this pipeline

````

---
### References
Under the hood this pipeline uses some coding tools, please include the following ciations in your work:

* Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. doi:10.1038/nbt.3820

* Team, R. C. (2017). R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna. http s. www. R-proje ct. org.

* Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686.

* Danecek, Petr, et al. "Twelve years of SAMtools and BCFtools." Gigascience 10.2 (2021): giab008.

* Suvakov, Milovan, et al. "CNVpytor: a tool for copy number variation detection and analysis from read depth and allele imbalance in whole-genome sequencing." Gigascience 10.11 (2021): giab074.  
---

### Contact
If you have questions, requests, or bugs to report, open an issue in github, or email <iaguilaror@gmail.com>

### Dev Team
Israel Aguilar-Ordonez <iaguilaror@gmail.com>   
Victor Trevino Alvarado <vtrevino@tec.mx>   
Eugenio Guzman Cerezo <eugenio.guzman@tec.mx>   

This code was developed as part of Israel Aguilar-Ordoñez’s postdoctoral research at Tecnológico de Monterrey during the 2024–2025. 

### Cite us
- TO-DO
