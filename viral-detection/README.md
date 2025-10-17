# viral hit detection
Repo with the code that was used to detect viral reads in oriGen

===  

- This pipeline is meant to show the code used in: TO-DO-add url and doi after paper is published

'viral-detection' is a pipeline tool that takes sequencing data in fastq format and a viral balst database to generate the following outputs:  

````
For multiple CRAM files:  

results/unmap_counts.tsv  # Table file with viral hits for all samples
results/dna_virus_analysis.ipynb => # Run this jupyter notebook to generate post-analysis plots and tables

````

---

### Features
  **-v 0.0.1**

* Supports CRAM files
* Results include TSV with viral hits, and many plots in a jupyter notebook
* Runs as bash scripts, and jupyter notebook

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
| [jupyterlab](https://anaconda.org/conda-forge/jupyterlab) | 4.4.9 | jupyter lab |
| [blast](https://anaconda.org/bioconda/blast) | 2.17.0 | makeblastdb, blastn |
| [parallel](https://anaconda.org/conda-forge/parallel) | 20250822 | parallel |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).  

#### python packages required:

```
numpy version: 1.26  
pandas version: 2.2.3  
```

---

### Installation
Download pipeline from Github repository:  
```
git clone https://github.com/Iaguilaror/oriGen-wgs001-scripts.git

cd oriGen-wgs001-scripts/viral-detection
```
---


## Replicate our analysis (Testing the pipeline):

* Estimated test time:  **10 minute(s)**  
* on a 16 core, 64G RAM machine  

1. To test pipeline execution using test data, run:  
```
bash runtest.sh
```

2. Your console should print the following message will appear:  
```
======
 Basic pipeline TEST SUCCESSFUL
======
```

3. Pipeline results for test data should be in the following directories:  
```
./results/
```
---


### Pipeline Inputs

* A directory containing many `CRAM files` with alignment data of many samples.  

Example contents for TWO samples  
```
test/
└── data
    ├── dummy2.cram
    ├── dummy2.cram.crai
    ├── dummy.cram
    └── dummy.cram.crai
```  

* A full set of `test/reference/virus_dumydb.*` files with the viral blastdb to survey.  

Example reference dir
```
test/reference/
├── virus_dumydb.fa
├── virus_dumydb.ndb
├── virus_dumydb.nhr
├── virus_dumydb.nin
├── virus_dumydb.njs
├── virus_dumydb.not
├── virus_dumydb.nsq
├── virus_dumydb.ntf
└── virus_dumydb.nto
```  

---

### Pipeline Results

```
Inside the directory sto/ you can find the following:  

results/
...
```


---

### module directory structure

````
.
.
├── 1_preprocess.sh             # first script to run. Pre-process CRAM files to blast results
├── dna_virus_analysis.ipynb    # Jupyter notebook for final processing
├── finds                       # temporary results
├── README.md                   # This README
├── scripts                     # bash and python script to process data
└── test                        # test data and references
````

---
### References
Under the hood this pipeline uses some coding tools, please include the following ciations in your work:

* Danecek, Petr, et al. "Twelve years of SAMtools and BCFtools." Gigascience 10.2 (2021): giab008.

---

### Contact
If you have questions, requests, or bugs to report, open an issue in github, or email <iaguilaror@gmail.com>

### Dev Team
Eugenio Guzman Cerezo <eugenio.guzman@tec.mx>   
Victor Trevino Alvarado <vtrevino@tec.mx>   
Israel Aguilar-Ordonez <iaguilaror@gmail.com>   

### Cite us
- TO-DO

