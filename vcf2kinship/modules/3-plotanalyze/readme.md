# plotanalyze  
**Author(s):**

* Israel Aguilar-Ordoñez (iaguilaror@gmail.com)

**Date:** July 2024  

---

## Module description:  

A (DSL2) Nextflow module to plot king2 results to identify related samples.


## Module Dependencies:
| Requirement | Version  | Required Commands |
|:---------:|:--------:|:-------------------:|
| [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 24.04.3 | nextflow |
| [R](https://www.r-project.org/) | 4.4.1 (2024-06-14) | Rscript |

#### R packages required:

```
cowplot version: 1.1.1
vroom version: 1.6.5
dplyr version: 1.1.4
ggplot2 version: 3.5.1
ggsci version: 3.2.0
patchwork version: 1.2.0
```

### Input(s):

* A `king.kin` file with king2 --related results.

* A `*_survey_completion.tsv*` file with a sample ID (same as original VCF) an a completion col (from 0 to 100).

Example contents  
```
==> test/data/king.kin <==
FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    HetConc HomIBS0 Kinship IBD1Seg IBD2Seg PropIBD InfType Error
ori     NA19648 NA19649 15715283        1.000   0.0000  0.0511  0.0284  0.2098  0.2187  -0.0194 0.0545  0.0000  0.0273  UN      0
ori     NA19648 NA19651 15713266        1.000   0.0000  0.0515  0.0280  0.2073  0.2168  -0.0151 0.0155  0.0000  0.0078  UN      0

==> test/reference/MXL_origen_dummy_survey_completion.tsv <==
ID      completion
NA19648 72
NA19649 67

```

### Outputs:

* Many `.png files` with plots showing relatedness results.  

Example contents  
```
1_allsamples.png                    # Shows all the inferred relations from PO(parent-offspring), FS(FullSiblings), 2nd, 3rd, 4th degree, and UNrelate.
2_related_samples.png               # same as 1, but only showing samples with at least 1(one) relation of any type.
3_types_of_pair.png                 # Summary of the type of relations found.
4_close_pairs.png                   # Same as 3, but only for PO, FS, or 2nd relations
5_samples_realtion_and_survey.png   # Number of relations by sample with PO,FS or at least 2nd. Also shows number of questionare completion.

```

## Module parameters:

NONE

## Testing the module:

* Estimated test time:  **1 minute(s)**  

1. Test this module locally by running,
```
bash testmodule.sh
```

2.`[>>>] Module Test Successful` should be printed in the console.  

## module directory structure

````
.
├── main.nf                                 # Nextflow main script
├── readme.md                               # this readme
├── scripts -> ../../scripts/               # symlink to general scripts directory
├── test                                    # dir with test materials
│   └── data -> ../../2-king/test/results/02-king/  # symlink to input data: .vcf.gz and .vcf.gz.tbi files go here
│   └── results                            
│       └── 03-plotking                                       
│           └── ...       # many plots
├── testmodule.nf                           # Nextflow test script to call the main.nf after simulating channel interactions
└── testmodule.sh                           # bash script to test the whole module
````
## References
* Team, R. C. (2017). R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna. http s. www. R-proje ct. org.

* Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686.