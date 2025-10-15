# oGpv3
Repo with the code that was used to call variants in oriGen

===  

- This pipeline is meant to show the code used in: TO-DO-add url and doi after paper is published

'oGpv3' is a pipeline tool that takes sequencing data in fastq format to generate the following outputs:  

````
For a set of 4 lanes, R1 and R2 pairs of fastq files

1) *.cram                     # Sample level CRAM file with alignments to GRCh38
2) *.vcf.gz (or bcf.gz)      # Sample level VCF file with SNV and INDEL variants
3) *.cnvpytor.bcf            # Sample level BCF file with copy number variants

````
---

### Features
  **-v 0.0.1**

* Supports FASTQ files (from illumina)
* Results include CRAM, VCF, and Copy Number BCF
* Scalability and reproducibility via a Nextflow-based framework   

---

## Requirements
#### This pipeline was successfully tested on the following OS*:  
* [Ubuntu 22.04.4 LTS](https://releases.ubuntu.com/focal/)

#### Incompatible OS*:
* UNKNOWN  

\* This pipeline may run in other UNIX based OS and versions, but testing is required.  

********* TODO: COMPLETE THIS *********

#### Command line Software required:
| Requirement | Version  | Required Commands * |
|:---------:|:--------:|:-------------------:|
| [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 24.04.3 | nextflow |
| [bcftools](https://anaconda.org/bioconda/bcftools) | 1.20 | bcftools |
| [vep](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) | 113.0 | vep |
| [R](https://www.r-project.org/) | 4.4.1 (2024-06-14) | Rscript |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).  

#### R packages required:

```
vroom version: 1.6.5
tidyr version: 1.3.1
dplyr version: 1.1.4
ggplot2 version: 3.5.1
purrr version: 1.0.4
ggalluvial version: 0.12.5
ggsci version: 3.2.0
scales version: 1.4.0
stringr version: 1.5.1
viridis version: 0.6.5
hexbin version: 1.28.5
ggExtra version: 0.10.1
patchwork version: 1.3.0
cowplots version: 1.1.3
```

---

### Installation
Download pipeline from Github repository:  
```
git clone https://github.com/Iaguilaror/oriGen-wgs001-scripts.git

cd oriGen-wgs001-scripts/vep-extended
```
---

### IMPORTANT: Install References

This pipeline requires you to download external databases and pass them to the runtest.sh script as the --dbsnp_ref --gnomad_ref --mcps_ref params.

It also requires that you have the following files/databases in your local system:

```
GRCh38_dbsnp156.vcf.gz
gnomAD_v4_for_VEP_annotation.vcf.gz
MCPS_for_VEP_annotation.vcf.gz 
```

Each database must have the same chr nomenclature as your VCF file.  

Db's were downloaded from: https://www.ncbi.nlm.nih.gov/snp/  , https://gnomad.broadinstitute.org/downloads  ,  https://www.ctsu.ox.ac.uk/research/prospective-blood-based-study-of-150-000-individuals-in-mexico    

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

* A directory containing a `.vcf.gz file and .vcf.gz.tbi index` with genotypes of multiple samples.

Example contents  
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
...
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA19648 NA19649 NA19650 NA19651 NA19652 NA19653 NA19654 NA19655 NA19656NA19657 NA19658 NA19659 NA19660 NA19661 NA19662 NA19663 NA19664 NA19665 NA19669 NA19670 NA19671 NA19675 NA19676 NA19677 NA19678 NA19679 NA19680NA19681 NA19682 NA19683 NA19684 NA19685 NA19686 NA19716 NA19717 NA19718 NA19719 NA19720 NA19721 NA19722 NA19723 NA19724 NA19725 NA19726 NA19727NA19728 NA19729 NA19730 NA19731 NA19732 NA19733 NA19734 NA19735 NA19740 NA19741 NA19746 NA19747 NA19748 NA19749 NA19750 NA19751 NA19752 NA19755NA19756 NA19757 NA19758 NA19759 NA19760 NA19761 NA19762 NA19763 NA19764 NA19770 NA19771 NA19772 NA19773 NA19774 NA19775 NA19776 NA19777 NA19778NA19779 NA19780 NA19781 NA19782 NA19783 NA19784 NA19785 NA19786 NA19787 NA19788 NA19789 NA19790 NA19792 NA19794 NA19795 NA19796
1       79279   .       C       T       .       PASS    AF=0.000343643;AN=172;AC=2;MAF=0.000333 GT      0/0     0/0     0/0     0/0     0/0
     ./.     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/1     0/1     ./.     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0
     0/0     ./.     ./.     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
1       86065   .       G       C       .       PASS    AF=0.041595;AN=174;AC=2;MAF=0.0403065   GT      0/0     0/0     0/0     0/0     0/0
     ./.     0/0     0/0     0/0     0/1     0/0     0/1     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0
     0/0     ./.     ./.     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
...
```  

* A `GRCh38_dbsnp156.vcf.gz` VCF file containing the dbSNP Build 156 database aligned to the human genome reference GRCh38.  

* A `gnomAD_v4_for_VEP_annotation.vcf.gz` VCF file containing the gnomAD v4 database (all chromosomes).  

* A `MCPS_for_VEP_annotation.vcf.gz` VCF file containing the Mexico City Prospective Project allele frequencies (all chromosomes).  

---

### Pipeline Results

```
Inside the directory test/results/origen-vepextended/05-0-extractnovel you can find the following:  

1) *.novel.vcf.gz                           # VEP annotated VCF file only for variants with no rsID found  

Inside the directory test/results/origen-vepextended/05-a-dbsnp_summ you can find the following:  

2) *.dbsnp_summary.txt                      # Summary of novel variants by type, and chromosome  

Inside the directory test/results/origen-vepextended/05-b-vcf2tsv you can find the following:  

3) *_gnomADv4_MCPS.allelefrequencies.tsv    # Table format for allele frecuencies in gnomAD, and MCPS  

Inside the directory test/results/origen-vepextended/06-a-dbsnp_plot you can find the following:  

4) novel_variants.tsv                       # Table with total % of novel variants by type  
5) sankey_novel.png                         # Plot for the novel_variants.tsv table  

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

* Karczewski, Konrad J., et al. "The mutational constraint spectrum quantified from variation in 141,456 humans." Nature 581.7809 (2020): 434-443.

* McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F.
The Ensembl Variant Effect Predictor. Genome Biology Jun 6;17(1):122. (2016)

* Ziyatdinov, Andrey, et al. "Genotyping, sequencing and analysis of 140,000 adults from Mexico City." Nature 622.7984 (2023): 784-793.  

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
