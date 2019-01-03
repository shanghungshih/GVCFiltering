# GVCFiltering
Filtering vcf(variant calling) with population gvcf(genome.vcf) from joint variant calling

## Prerequisite
- JDK 11

## Installation
``` shell
git clone https://github.com/shanghungshih/GVCFiltering.git
```

## Quick start
### Filtering with Allele Frequency(AF) and Confidence(Conf)
- input 4 arguments:
    - AF: keep variants of vcf which AF <= $AF (ex. `0.01`, which means to keep "the variants with at least 1% of population appearance")
    - Conf: keep variants filtered by MAF which Conf > $Conf (ex. `0.8`, which means to keep "the variants with genotype information > 80% of population")
        - Calculation of confidence for each variant in population gvcf: # of samples with genotype information at the loci / # of samples
    - Population gvcf file: (ex. `population_test.gvcf`)
    - vcf file to be filter: (ex. `test.vcf`) 
```
java -jar GVCFiltering.jar 0.01 0.8 population_test.gvcf test.vcf
```

