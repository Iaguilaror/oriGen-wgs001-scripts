# origen-LoFbySample
Repo to annotate LoF variants, and count them by sample

===  

- A tool for running VEP LoF (Loss of Function annotation) pluggin in VCF files

- This pipeline is meant to reproduce the results in: TO-DO-add url and doi after paper is published

'LoFbySample' is a pipeline tool that takes variant data in VCF format, and loftee databases to to generate the following outputs:  

````
For ONE VCF and .tbi pair of files

1) *_clean_vep_loftee_only.tsv          # variants annotated with VEP, including  LoF     LoF_filter      LoF_flags       LoF_info fields   
2) *_clean_vep_loftee_only_summary.tsv  # a summary by sample of homozygous and heterozygous LoF variants
3) n_HC_LoF_variants.svg                # plot showing the distribution of LoF variants in all samples
4) summary_stats_LoF_allsamples.tsv     # a table with mean, sd, se, min, max, etc number of LoF vartiants by type

````
---

### Features
  **-v 0.0.1**

* Supports VCF and .tbi files ##fileformat=VCFv4.2
* Results include plots that show distribution of LoF variation
* Results include tables with LoF annotation
* Scalability and reproducibility via a Nextflow-based framework   

---

## Requirements
#### Compatible OS*:
* [Ubuntu 22.04.4 LTS](https://releases.ubuntu.com/focal/)

#### Incompatible OS*:
* UNKNOWN  

\* origen-LoFbySample may run in other UNIX based OS and versions, but testing is required.  

#### Command line Software required:
| Requirement | Version  | Required Commands * |
|:---------:|:--------:|:-------------------:|
| [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 24.04.3 | nextflow |
| [bcftools](https://anaconda.org/bioconda/bcftools) | 1.20 | bcftools |
| [vep](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) | 113.0 | vep |
| [LoF vep plugin](https://github.com/konradjk/loftee) | 1.0.4 | NA |
| [R](https://www.r-project.org/) | 4.4.1 (2024-06-14) | Rscript |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).  

#### R packages required:

```
vroom version: 1.6.5
tidyr version: 1.3.1
dplyr version: 1.1.4
ggplot2 version: 3.5.1
```

---

### Installation
Download pipeline from Github repository:  
```
git clone https://github.com/Iaguilaror/oriGen-wgs001-scripts.git

cd oriGen-wgs001-scripts/LoFbySample
```
---

### IMPORTANT: Install References

This pipeline requires you to install the loftee VEP plugin following this instructions: https://github.com/konradjk/loftee    

It also requires that you have the following files/databases in your local system:

```
human_ancestor_GRCh38_autosomes.fa.gz
loftee.sql
```
Follow the https://github.com/konradjk/loftee instructions to download them.  

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

* A `human_ancestor_GRCh38_autosomes.fa.gz` FASTA-format file containing DNA sequences of inferred ancestral human genome regions corresponding to the autosomal chromosomes (1–22) of the GRCh38 reference assembly. See https://github.com/konradjk/loftee  

* A `loftee.sql` structured query language (SQL) script that defines or populates database tables used internally by LOFTEE to support its annotation logic. See https://github.com/konradjk/loftee   

---

### Pipeline Results

```
Inside the directory test/results/origen-LoFbysample-results/03-filterlof you can find the following:

1) *_clean_vep_loftee_only.tsv          # variants annotated with VEP, including  LoF     LoF_filter      LoF_flags       LoF_info fields   

Inside the directory test/results/origen-LoFbysample-results/04-countlof you can find the following:

2) *_clean_vep_loftee_only_summary.tsv  # a summary by sample of homozygous and heterozygous LoF variants

Inside the directory test/results/origen-LoFbysample-results/05-gather_plot you can find the following:

3) n_HC_LoF_variants.svg                # plot showing the distribution of LoF variants in all samples
4) summary_stats_LoF_allsamples.tsv     # a table with mean, sd, se, min, max, etc number of LoF vartiants by type

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
