# origen-VCF2PCP-2
NF Pipeline to perform PCA and ADMIXTURE in origen data

===  

- A tool for running PCA and ADMIXTURE in VCF files

- This pipeline is meant to reproduce the results in: TO-DO-add url and doi after paper is published

'vcf2pcp' is a pipeline tool that takes variant data in VCF format, and sample population annotation to generate the followin outputs:  

````
For ONE VCF and .tbi pair of files, and the corresponding sample population annotation TSV

1) many bipca.PCXvsPCY.svg          # for example, PC1vsPC2 or PC1vsPC3. A 2-axis PCA plot.
2) parallel_plot.svg                # Parallel Coordinate Plot summarizing all the samples, by annotated population
3) parallel_plot.PCA_df.tsv         # PCA coordinates for the first 10 PComponents, in table format
4) parallel_plot.significant_pc.tsv # explained_variance_proportions for each PC
5) allwgs.converted2plink.5.admixture_plot.svg # SVG figure with column plot for admixture proportions. One for each K calculated
6) allwgs.converted2plink.5.admixture_proportion.tsv # Table to create the plot in allwgs.converted2plink.5.admixture_plot.svg. One for each K calculated
7) all_admixture.svg                # Big figure with all the admixture plots for all K calculated.

````
---

### Features
  **-v 0.0.1**

* Supports VCF and .tbi files ##fileformat=VCFv4.2
* Results include plots that show proportions of amixture
* Calculates only the first 10 PC
* Results include plots that show many PCA dimensions
* Results include tables to build the ADMIXTURE and PCA plots
* Scalability and reproducibility via a Nextflow-based framework   

---

## Requirements
#### Compatible OS*:
* [Ubuntu 22.04.4 LTS](https://releases.ubuntu.com/focal/)

#### Incompatible OS*:
* UNKNOWN  

\* origen-vcf2kinship may run in other UNIX based OS and versions, but testing is required.  

#### Command line Software required:
| Requirement | Version  | Required Commands * |
|:---------:|:--------:|:-------------------:|
| [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 24.04.3 | nextflow |
| [plink2](https://anaconda.org/bioconda/plink2) | v2.00a5.12LM | plink2 |
| [eigensoft](https://anaconda.org/bioconda/eigensoft) | 8.0.0 | smartpca |
| [admixture](https://anaconda.org/bioconda/admixture) | 1.3.0 | admixture |
| [R](https://www.r-project.org/) | 4.4.1 (2024-06-14) | Rscript |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).  

#### R packages required:

```
vroom version: 1.6.5
tidyr version: 1.3.1
dplyr version: 1.1.4
ggplot2 version: 3.5.1
ggsci version: 3.2.0
patchwork version: 1.2.0
svglite version: 2.2.0
cowplot version: 1.1.3
scales version: 1.4.0
stringr version: 1.5.1
RColorBrewer version: 1.1.3

```


---

### Installation
Download pipeline from Github repository:  
```
git clone https://github.com/Iaguilaror/oriGen-wgs001-scripts.git

cd oriGen-wgs001-scripts/vcf2pcp
```

---

## Replicate our analysis (Testing the pipeline):

* Estimated test time:  **10 minute(s)**  

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

* A directory containing a `.vcf.gz file and .vcf.gz.tbi index` with genotypes of multiple samples.

Example contents  
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
...
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096 HG00097 HG00099 HG00100
21      5225159 .       G       C       .       .       AN=5000;AF=0.0382141;MAF=0.0382141;AC=190       GT      0/0     0/0     0/0     0/0
21      5244487 .       A       C       .       .       AN=5000;AF=0.0211014;MAF=0.0211014;AC=115       GT      0/0     0/0     0/0     0/0
21      5254932 .       G       A       .       .       AN=5000;AF=0.771873;MAF=0.228127;AC=3561        GT      0/1     1/1     1/1     1/1
21      7325410 .       C       T       .       .       AN=5000;AF=0.0183994;MAF=0.0183994;AC=135       GT      0/0     0/0     0/0     0/0
...
```  

* A `sample_annotations_test.tsv` file with a sample ID (same as original VCF) an population annotations.

Example contents  
```
sample  Population code Superpopulation code
HG00096 GBR     EUR
HG00097 ORI     ORI
HG00099 ORI     ORI
HG00100 GBR     EUR
...
```  

* A `refpops-for-projection.txt` file with Superpopulation code's (same as sample_annotations_test.tsv). All samples from these populations will be used as 'references' for ADMIXTURE and PCA. Any other sample will be 'projected' on the data from these refpops.

Example contents  
```
AFR
EUR
SAS
EAS
AMR
...
```  
---

### Pipeline Results

```
Inside the directory test/results/origen-VCF2PCP-2-results/04b-projected-plotpca you can find the following:

1) many bipca.PCXvsPCY.svg          # for example, PC1vsPC2 or PC1vsPC3. A 2-axis PCA plot.
2) parallel_plot.svg                # Parallel Coordinate Plot summarizing all the samples, by annotated population
3) parallel_plot.PCA_df.tsv         # PCA coordinates for the first 10 PComponents, in table format
4) parallel_plot.significant_pc.tsv # explained_variance_proportions for each PC

Inside the directory test/results/origen-VCF2PCP-2-results/05b-plot_admixture/ you can find the following:

5) allwgs.converted2plink.5.admixture_plot.svg # SVG figure with column plot for admixture proportions. One for each K calculated
6) allwgs.converted2plink.5.admixture_proportion.tsv # Table to create the plot in allwgs.converted2plink.5.admixture_plot.svg. One for each K calculated

Inside the directory test/results/origen-VCF2PCP-2-results/05c-gather_admixture/ you can find the following:

7) all_admixture.svg                # Big figure with all the admixture plots for all K calculated.

```


---

### module directory structure

````
.
├── configfiles         # the Nextflow configuration files to tune runs for different hardware
├── main.nf         # the Nextflow main script
├── modules/        # sub-dirs for development of the Nextflow modules
├── README.md       # This readme
├── runtest.sh      # bash script to launch the pipeline test locally
├── scripts/        # directory with all the scripts used by the pipeline
└── test
    └── data       # sample data to run this pipeline
    └── reference  # sample annotation or reference files to test this pipeline

````

---
### References
Under the hood Proteomic compare uses some coding tools, please include the following ciations in your work:

* Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. doi:10.1038/nbt.3820

* Team, R. C. (2017). R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna. http s. www. R-proje ct. org.

* Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686.

* Danecek, Petr, et al. "Twelve years of SAMtools and BCFtools." Gigascience 10.2 (2021): giab008.

* Purcell, Shaun, et al. "PLINK: a tool set for whole-genome association and population-based linkage analyses." The American journal of human genetics 81.3 (2007): 559-575.

* Manichaikul, Ani, et al. "Robust relationship inference in genome-wide association studies." Bioinformatics 26.22 (2010): 2867-2873.


---

### Contact
If you have questions, requests, or bugs to report, open an issue in github, or email <iaguilaror@gmail.com>

### Dev Team
Israel Aguilar-Ordonez <iaguilaror@gmail.com>   

### Cite us
- TO-DO
