# vcf2plink  
**Author(s):**

* Israel Aguilar-Ordoñez (iaguilaror@gmail.com)

**Date:** July 2024  

---

## Module description:  

A (DSL2) Nextflow module to convert vcf files into plink format (bim bed fam).


## Module Dependencies:
| Requirement | Version  | Required Commands |
|:---------:|:--------:|:-------------------:|
| [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 24.04.3 | nextflow |
| [plink2](https://anaconda.org/bioconda/plink2) | v2.00a5.12LM | plink2 |

### Input(s):

* A `.vcf.gz file and .vcf.gz.tbi index` with genotypes of multiple samples.

Example contents  
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
...
CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA19648 NA19649 NA19651 NA19652       NA19654 NA19655 NA19657 NA19658 NA19660 NA19661 NA19663 NA19664 NA19669 NA19670 NA19676 NA19678       NA19679 NA19681 NA19682 NA19684 NA19716 NA19717 NA19719 NA19720 NA19722 NA19723 NA19725 NA19726       NA19728 NA19729 NA19731 NA19732 NA19734 NA19735 NA19740 NA19741 NA19746 NA19747 NA19749 NA19750       NA19752 NA19755 NA19756 NA19758 NA19759 NA19761 NA19762 NA19764 NA19770 NA19771 NA19773 NA19774       NA19776 NA19777 NA19779 NA19780 NA19782 NA19783 NA19785 NA19786 NA19788 NA19789 NA19792 NA19794       NA19795
chr1    10150   .       C       T       .       PASS    AF=0.0109375;AN=14;AC=1;MAF=0.0145278   GT   ./.      ./.     ./.     ./.     ./.     0/0     0/0     0/0     ./.     ./.     ./.     ./.     ./.  ./.      ./.     ./.     ./.     ./.     0/0     ./.     ./.     ./.     ./.     ./.     ./.     ./.  ./.      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     0/0     ./.     ./.     0/0  0/1      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.  ./.      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.
chr1    10327   .       T       C       .       PASS    AF=0.250478;AN=34;AC=9;MAF=0.212662     GT   0/0      ./.     0/1     1/1     ./.     ./.     ./.     0/1     ./.     ./.     0/0     ./.     ./.  ./.      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     1/1     ./.  ./.      ./.     0/1     ./.     0/0     ./.     ./.     0/1     ./.     0/0     ./.     ./.     ./.  0/0      0/1     ./.     0/0     ./.     ./.     ./.     ./.     0/0     ./.     0/0     ./.     0/0  ./.      ./.     ./.     ./.     ./.     ./.     ./.     ./.     0/0     ./.     ./.     ./.     ./.
...
```

### Outputs:

* The `.converted2plink.bed .converted2plink.bim .converted2plink.fam files` with GT in plink format.  

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
│   └── data -> ../../0-keepautosomes/test/results/0-keepautosomes/ # symlink to input data: .vcf.gz and .vcf.gz.tbi files go here
│   └── results                            
│       └── 01-vcf2plink                                       
│           └── *.converted2plink.bed       # variants in plink format
│           └── *.converted2plink.bim       # must come with the .bed file
│           └── *.converted2plink.fam       # must come with the .bed file
├── testmodule.nf                           # Nextflow test script to call the main.nf after simulating channel interactions
└── testmodule.sh                           # bash script to test the whole module
````
## References
* Purcell, Shaun, et al. "PLINK: a tool set for whole-genome association and population-based linkage analyses." The American journal of human genetics 81.3 (2007): 559-575.