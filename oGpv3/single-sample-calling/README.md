# oGpv3
Repo with the code that was used to call variants in oriGen

===  

- This pipeline is meant to show the code used in: TO-DO-add url and doi after paper is published

'oGpv3' is a pipeline tool that takes sequencing data in fastq format to generate the following outputs:  

````
For a set of 4 lanes, R1 and R2 pairs of fastq files:  

sto/BCFS/CNV/MYRUNID1/MYCODEID1.bcf  # Sample level BCF file with copy number variants  
sto/BCFS/SNV/MYRUNID1/MYCODEID1.bcf  # Sample level BCF file with SNV and INDEL variants  
sto/CRAMS/MYRUNID1/MYCODEID1.cram    # Sample level CRAM file with alignments to GRCh38  
````

---

### Features
  **-v 0.0.1**

* Supports FASTQ files (from illumina)
* Results include CRAM, short variant BCF, and Copy Number BCF
* Scalability and reproducibility via a snakemake-based framework   

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
| [snakemake](https://anaconda.org/bioconda/snakemake) | 9.13.2 | snakemake |
| [fastp](https://anaconda.org/bioconda/fastp) | 1.0.1 | fastp |
| [bwa-mem2](https://anaconda.org/bioconda/bwa-mem2) | 2.3 | bwa-mem2 |
| [samtools](https://anaconda.org/bioconda/samtools) | 1.22.1 | samtools |
| [graphtyper](https://anaconda.org/bioconda/graphtyper) | 2.7.7 | graphtyper |
| [cnvpytor](https://anaconda.org/bioconda/cnvpytor) | 1.3.1 | cnvpytor |
| [python](https://anaconda.org/anaconda/python/files?version=3.11.4) | 3.11 | python |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).  

#### python packages required:

```
numpy version: 1.26  
```

---

### Installation
Download pipeline from Github repository:  
```
git clone https://github.com/Iaguilaror/oriGen-wgs001-scripts.git

cd oriGen-wgs001-scripts/oGpv3/single-sample-calling
```
---


## Replicate our analysis (Testing the pipeline):

* Estimated test time:  **10 minute(s)**  
* on a 16 core, 64G RAM machine  

1. To test pipeline execution using test data, run:  
```
bash runtest.sh
```

2. Your console should print the snakemake log for the run, once every process has been submitted, the following message will appear:  
```
======
 Basic pipeline TEST SUCCESSFUL
======
```

3. Pipeline results for test data should be in the following directories:  
```
./sto/BCFS
./sto/CRAMS
```
---


### Pipeline Inputs

* A directory containing a set of `multi LANE fastq files` with sequencing data of one sample. However, many samples can be run in parallel

Example contents for ONE sample  
```
input/
└── MYRUNID1
   └── MYCODEID1
       ├── MYCODEID1_A0_L001_R1_000.fastq.gz
       ├── MYCODEID1_A0_L001_R2_000.fastq.gz
       ├── MYCODEID1_A0_L002_R1_000.fastq.gz
       ├── MYCODEID1_A0_L002_R2_000.fastq.gz
       ├── MYCODEID1_A0_L003_R1_000.fastq.gz
       ├── MYCODEID1_A0_L003_R2_000.fastq.gz
       ├── MYCODEID1_A0_L004_R1_000.fastq.gz
       └── MYCODEID1_A0_L004_R2_000.fastq.gz
```  

* A `ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa` fasta file cotaining the human genome reference GRCh38.  

* A full set of bwa-mem2 references created from the above fasta file.  

Example reference dir
```
ref/
├── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
├── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.0123
├── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.amb
├── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.ann
├── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.bwt.2bit.64
├── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai
└── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.pac
```  

---

### Pipeline Results

```
Inside the directory sto/ you can find the following:  

sto/
├── BCFS
│   ├── CNV
│   │   └── MYRUNID1
│   │       ├── MYCODEID1.bcf      # Copy Number Variants detected by CNVpytor
│   │       └── MYCODEID1.bcf.csi
│   └── SNV
│       └── MYRUNID1
│           ├── MYCODEID1.bcf      # Short Variants (SNV, and INDL) detecetd by graphtyper
│           └── MYCODEID1.bcf.csi
└── CRAMS
    └── MYRUNID1
        ├── MYCODEID1.cram         # Alignment files
        └── MYCODEID1.cram.crai
```


---

### module directory structure

````
.
├── anno         # not used
├── blob         # not used
├── cache        # not used
├── input        # directory with FASTQ files in subdirectories by sample
├── logs         # snakemake logs
├── oGpv3.smk    # the snakemake main script
├── README.md    # This readme
├── ref          # directory to store the human reference files (fasta, and bwa-mem2 index files)
├── reports      # not used
├── runtest.sh   # script to run the test data flow
└── sto          # storage directory for results
````

---
### References
Under the hood this pipeline uses some coding tools, please include the following ciations in your work:

* Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.  

* Danecek, Petr, et al. "Twelve years of SAMtools and BCFtools." Gigascience 10.2 (2021): giab008.

* Eggertsson, Hannes P., et al. "Graphtyper enables population-scale genotyping using pangenome graphs." Nature genetics 49.11 (2017): 1654-1660.  

* Suvakov, Milovan, et al. "CNVpytor: a tool for copy number variation detection and analysis from read depth and allele imbalance in whole-genome sequencing." Gigascience 10.11 (2021): giab074.  
---

### Contact
If you have questions, requests, or bugs to report, open an issue in github, or email <iaguilaror@gmail.com>

### Dev Team
Eugenio Guzman Cerezo <eugenio.guzman@tec.mx>   
Victor Trevino Alvarado <vtrevino@tec.mx>   
Israel Aguilar-Ordonez <iaguilaror@gmail.com>   

### Cite us
- TO-DO
