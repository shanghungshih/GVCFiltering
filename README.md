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
- required 4 arguments:
    - `--allele_frequency` `-af`: Keep variants with allele frequency <= threshold (ex. `0.01`, which means to keep "the variants with at least 1% of population appearance")
    - `--confidence` `-conf`: Keep variants with confidence > threshold (ex. `0.8`, which means to keep "the variants with genotype information > 80% of population")
        - Calculation of confidence for each variant in population gvcf: # of samples with genotype information at the loci / # of samples
    - `--population` `-p`: Input population gvcf file (ex. `population_test.gvcf`)
    - `--vcf` `-v`: Input vcf file to be filtered (ex. `test.vcf`)

- optional argument:
    - `--export_population` `-e`: Export population gvcf file with GTF(genotype frequency), AF(allele frequency) and CONF(confidence) annotation
    
```
java -jar GVCFiltering.jar -af 0.01 -conf 0.8 -p population_test.gvcf -v test.vcf -e
```

### output
- output.vcf (Note: one alternative per line!)
- output.log
- (option) output.g.vcf
