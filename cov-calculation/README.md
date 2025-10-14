# origen-cov-calculation
Small pipeline to calculate coverage in genes in the oriGen dataset  

===  

- A tool for analyzing CRAM files for coverage  

- This pipeline is meant to reproduce the results in: TO-DO-add url and doi after paper is published  

'cov-calculation' is a pipeline tool that takes CRAM alignment data to generate the following outputs:  

````
For MULTIPLE CRAM and .crai pair of files

1) samtools_mean_dp_all_cov_for_cnv.tsv.gz  # Table format with mean cov for each exon, for all samples 

````
---

### Features
  **-v 0.0.1**

* Supports CRAM and .crai files
* Results include a table with All exon coverage by sample, all samples
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
| [bedtools](https://anaconda.org/bioconda/bedtools) | 2.31.1 | bedtools |
| [samtools](https://anaconda.org/bioconda/samtools) | 1.22.1 | samtools |
| [R](https://www.r-project.org/) | 4.4.1 (2024-06-14) | Rscript |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).  

#### R packages required:

```
-vroom version: 1.6.5
-tidyr version: 1.3.1
-dplyr version: 1.1.4
-purr version: 1.0.4
-stringr version: 1.5.1
-biomaRt version: 2.62.1
```

---

### Installation
Download pipeline from Github repository:  
```
git clone https://github.com/Iaguilaror/oriGen-wgs001-scripts.git

cd oriGen-wgs001-scripts/cov-calculation
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

* A directory containing many `.cram file and .crai index` with GRCh38 alignment for many samples

Example contents  
```
HD     VN:1.6  SO:coordinate
@SQ     SN:1    LN:248956422    M5:2648ae1bacce4ec4b6cf337dcae37816     UR:/home/bioinfouser/references/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
A01314:63:H7VJTDSX7:2:1506:29948:8766   99      1       103549894       60      150M    =       103549974       230     GGTTCATGCCTGTAATCCTAGCACATTGGGAGGCCGAAGCGGGGGGGATCACTTGAGGTCAGGAGTTCAAGACCAGCCTGGCCTATATAGTGAAACCCCATCTCTACTAAAAATACAAAAATTAGCCAAGCATGGCGGTAGACGGCTGTA   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF   AS:i:150        XS:i:56 MQ:i:60 MC:Z:150M        ms:i:5514       MD:Z:150        NM:i:0  RG:Z:030L2
A01314:63:H7VJTDSX7:2:1155:5466:23140   83      1       103549910       60      149M    =       103549683       -376    CCTAGCACATTGGGAGGCCGAAGCGGGGGGGATCACTTGAGGTCAGGAGTTCAAGACCAGCCTGGCCTATATAGTGAAACCCCATCTCTACTAAAAATACAAAAATTAGCCAAGCATGGCGGTAGACGGCTGTAATCCCAGCTACCCGG    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    AS:i:149        XS:i:56 MQ:i:60 MC:Z:151M        ms:i:5491       MD:Z:149        NM:i:0  RG:Z:030L2
...
```  

---

### Pipeline Results

```
Inside the directory test/results/origen-cov-calculation-results/02-gatherchunks you can find the following:  


1) samtools_mean_dp_all_cov_for_cnv.tsv.gz  # Table format with mean cov for each exon, for all samples 

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

