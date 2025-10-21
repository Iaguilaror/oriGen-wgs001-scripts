# origen-vcf2kinship
NF pipeline to estimate kinship coefficients with King2 for all pairwise relationships in the origen data

===  

- A tool for finding related samples in VCF files

- This pipeline is meant to reproduce the results in: TO-DO-add url and doi after paper is published

The pipeline takes as INPUT a VCF and .tbi file created by joint genotyping, with ONLY SNPs (No indels). The end result is a set of plots indicating a set of related samples. It also gives you a VCF file with unrelated samples up to 2nd degree.

'vcf2kin' is a pipeline tool that takes variant data in VCF format and process it to generate the followin outputs:  

````
For ONE VCF and .tbi pair of files, you get these outputs

1) 1_allsamples.png                    # Shows all the inferred relations from PO(parent-offspring), FS(FullSiblings), 2nd, 3rd, 4th degree, and UNrelate.
2) 2_related_samples.png               # same as 1, but only showing samples with at least 1(one) relation of any type.
3) 3_types_of_pair.png                 # Summary of the type of relations found.
4) 4_close_pairs.png                   # Same as 3, but only for PO, FS, or 2nd relations
5) 5_samples_realation.png             # Number of relations by sample with PO,FS or at least 2nd. 
6) 6_IBD_ternaryplot.png               # Triangular plot for IBD0 IBD1 IBD2 coordinates
7) 7_Distribution_by-number_of_relatives.png   # Columns plot to show how many participants have N relatives
8) 8a_relatedness_bar.png              # Column plot showing % of participants related up to 3rd degree
9) 8b_close_relatedness_bar.png        # Column plot showing % of participants related up to 2nd degree
10) king.kin                           # KING2 results file. Shows pairwise IBD calcs and relationship inferred type 
11) barra_related.png                  # Column plot showing N and % of samples that will be removed.
12) samples_to_remove.txt              # list of ID samples that will be removed to eliminate kin relationships up to 2nd degree
13) *.only_unrelated_samples.vcf.gz    # VCF file, without samples related up to 2nd degree
````

---

### Features
  **-v 0.0.1**

* Supports VCF and .tbi files ##fileformat=VCFv4.2
* Results include plots that identify the related samples
* Results include a VCF without related samples up to 2nd degree
* Scalability and reproducibility via a Nextflow-based framework   

---

## Requirements
#### This pipeline was successfully tested on the following OS*:  
* [Ubuntu 22.04.4 LTS](https://releases.ubuntu.com/focal/)

#### Incompatible OS*:
* UNKNOWN  

\* This pipeline may run in other UNIX based OS and versions, but testing is required.  

#### Command line Software required:
| Requirement | Version  | Required Commands * |
|:---------:|:--------:|:-------------------:|
| [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 24.04.3 | nextflow |
| [bcftools](https://anaconda.org/bioconda/bcftools) | 1.20 | bcftools |
| [plink2](https://anaconda.org/bioconda/plink2) | v2.00a5.12LM | plink2 |
| [King](https://anaconda.org/bioconda/king) | 2.2.7 | king |
| [R](https://www.r-project.org/) | 4.4.1 (2024-06-14) | Rscript |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).  

#### R packages required:

```
vroom version: 1.6.5
dplyr version: 1.1.4
ggplot2 version: 3.5.1
ggsci version: 3.2.0
patchwork version: 1.2.0
tidyverse version: 2.0.0
scales version: 1.4.0
ggtern version: 3.5.0
```

---

### Installation
Download pipeline from Github repository:  
```
git clone https://github.com/Iaguilaror/oriGen-wgs001-scripts.git

cd oriGen-wgs001-scripts/vcf2kinship
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

* A directory containing a `.vcf.gz file and .vcf.gz.tbi index` with genotypes of multiple samples.

* VCF file must have 1,2,3.. chromosome name format (NOT chr1,chr2,chr3...)

Example contents  
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
...
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA19648 NA19649 NA19650 NA19651 NA19
652     NA19653 NA19654 NA19655 NA19656 NA19657 NA19658 NA19659 NA19660 NA19661 NA19662 NA19663 NA19664 NA19
665     NA19669 NA19670 NA19671 NA19675 NA19676 NA19677 NA19678 NA19679 NA19680 NA19681 NA19682 NA19683 NA19
684     NA19685 NA19686 NA19716 NA19717 NA19718 NA19719 NA19720 NA19721 NA19722 NA19723 NA19724 NA19725 NA19
726     NA19727 NA19728 NA19729 NA19730 NA19731 NA19732 NA19733 NA19734 NA19735 NA19740 NA19741 NA19746 NA19
747     NA19748 NA19749 NA19750 NA19751 NA19752 NA19755 NA19756 NA19757 NA19758 NA19759 NA19760 NA19761 NA19
762     NA19763 NA19764 NA19770 NA19771 NA19772 NA19773 NA19774 NA19775 NA19776 NA19777 NA19778 NA19779 NA19
780     NA19781 NA19782 NA19783 NA19784 NA19785 NA19786 NA19787 NA19788 NA19789 NA19790 NA19792 NA19794 NA19
795     NA19796
1       79279   .       C       T       .       PASS    AF=0.000343643;AN=172;AC=2;MAF=0.000333 GT      0/0
        0/0     0/0     0/0     0/0     ./.     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0
        0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
        0/0     0/0     0/0     ./.     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     ./.
        0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0         0/0     0/0     0/0     0/0     0/0     0/0     0/1     0/1     ./.     0/0     0/0     0/0     0/0         0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     ./.     ./.         0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0         0/0     0/0     0/0     0/0     0/0
1       86065   .       G       C       .       PASS    AF=0.041595;AN=174;AC=2;MAF=0.0403065   GT      0/0         0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/1     0/0     0/1     0/0     0/0         0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0         0/0     0/0     0/0     ./.     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     ./.         0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0         0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0         0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     ./.     ./.         0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0         0/0     0/0     0/0     0/0     0/0
...
```  

---

### Pipeline Results

```
Inside the directory test/results/origen-kinship-results/03-plotking you can find the following:

1) 1_allsamples.png                    # Shows all the inferred relations from PO(parent-offspring), FS(FullSiblings), 2nd, 3rd, 4th degree, and UNrelate.
2) 2_related_samples.png               # same as 1, but only showing samples with at least 1(one) relation of any type.
3) 3_types_of_pair.png                 # Summary of the type of relations found.
4) 4_close_pairs.png                   # Same as 3, but only for PO, FS, or 2nd relations
5) 5_samples_realtion_and_survey.png   # Number of relations by sample with PO,FS or at least 2nd. Also shows number of questionare completion.   
6) 6_IBD_ternaryplot.png               # Triangular plot for IBD0 IBD1 IBD2 coordinates
7) 7_Distribution_by-number_of_relatives.png   # Columns plot to show how many participants have N relatives
8) 8a_relatedness_bar.png              # Column plot showing % of participants related up to 3rd degree
9) 8b_close_relatedness_bar.png        # Column plot showing % of participants related up to 2nd degree

Inside the directory test/results/03b-select2remove you can find the following:

11) barra_related.png                  # Column plot showing N and % of samples that will be removed.
12) samples_to_remove.txt              # list of ID samples that will be removed to eliminate kin relationships up to 2nd degree

Inside the directory test/results/04-getsamples-from-vcf you can find the following:

13) *.only_unrelated_samples.vcf.gz    # VCF file, without samples related up to 2nd degree
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

* Purcell, Shaun, et al. "PLINK: a tool set for whole-genome association and population-based linkage analyses." The American journal of human genetics 81.3 (2007): 559-575.

* D.H. Alexander, J. Novembre, and K. Lange. Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655–1664, 2009.

* Ma, Jianzhong, and Christopher I. Amos. "Theoretical formulation of principal components analysis to detect and correct for population stratification." PloS one 5.9 (2010): e12510.

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
