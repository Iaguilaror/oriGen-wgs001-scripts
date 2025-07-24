# king  
**Author(s):**

* Israel Aguilar-Ordoñez (iaguilaror@gmail.com)

**Date:** July 2024  

---

## Module description:  

A (DSL2) Nextflow module to Infer Relatedness by IBD and kinship, using King2.


## Module Dependencies:
| Requirement | Version  | Required Commands |
|:---------:|:--------:|:-------------------:|
| [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 24.04.3 | nextflow |
| [King2](https://anaconda.org/bioconda/king) | | king |

### Input(s):

* A set of `.bed .bim .fam files` with genotypes of multiple samples.

Example contents  
```
*.bed
BINARY FILE

*.bim
1       .       0       10150   T       C
1       .       0       10327   C       T
1       .       0       10330   A       C
1       .       0       10352   A       T

*.fam
ori     NA19648 0       0       0       -9
ori     NA19649 0       0       0       -9
ori     NA19650 0       0       0       -9
ori     NA19651 0       0       0       -9
```

### Outputs:

* The `kingallsegs.txt and  king.kin files` with relatedness results.  

Example contents  
```
==> test/results/02-king/kingallsegs.txt <==
Segment Chr     StartMB StopMB  Length  N_SNP   StartSNP        StopSNP
1       1       0.010   125.185 125.174 751910  .       .
2       1       143.185 248.946 105.761 573456  .       .

==> test/results/02-king/king.kin <==
FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    HetConc HomIBS0 Kinship IBD1Seg IBD2Seg PropIBD InfType Error
ori     NA19648 NA19649 15715283        1.000   0.0000  0.0511  0.0284  0.2098  0.2187  -0.0194 0.0545  0.0000  0.0273  UN      0
ori     NA19648 NA19651 15713266        1.000   0.0000  0.0515  0.0280  0.2073  0.2168  -0.0151 0.0155  0.0000  0.0078  UN      0
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
├── test                                    # dir with test materials
│   └── data -> ../../1-vcf2plink/test/results/01-vcf2plink/    # symlink to input data: .vcf.gz and .vcf.gz.tbi files go here
│   └── results                            
│       └── 02-king                                       
│           └── kingallsegs.txt       # all segments analyzed by king
│           └── king.kin              # per pair of samples results (see https://www.kingrelatedness.com/manual.shtml)
├── testmodule.nf                           # Nextflow test script to call the main.nf after simulating channel interactions
└── testmodule.sh                           # bash script to test the whole module
````
## References
* Manichaikul, Ani, et al. "Robust relationship inference in genome-wide association studies." Bioinformatics 26.22 (2010): 2867-2873.