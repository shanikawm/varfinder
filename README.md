# VarFinder
## Comprehensive Tool for Simulating Reads and Comparing the Read Mappings and Variant Calling Algorithms

## Installation
Make sure the dependencies are satisfied. Then copy the file released script file set. 
### Dependencies
**OS:** Any Unix version which satisfies the below dependencies. Tested with Enterprise Linux 7

1. [htslib, tabix, bgzip](https://github.com/samtools/htslib) (Tested with htslib 1.17)
2. [SAMtools](https://github.com/samtools/samtools) (Tested with Version: 1.17 (using htslib 1.17))
3. [BCFtools](https://github.com/samtools/bcftools) (Tested with Version: 1.17 (using htslib 1.17))
4. [wgsim](https://github.com/lh3/wgsim) (Tested with Version: 1.17)
5. [ngsngs](https://github.com/RAHenriksen/NGSNGS) (Tested with version 0.9.0)
6. [bwa](https://github.com/lh3/bwa) (Tested with version 0.7.17-r1198-dirty)
7. [bowtie2](https://github.com/BenLangmead/bowtie2) (Tested with version 2.5.1)
8. [vg](https://github.com/vgteam/vg) (Tested with version v1.51.0 "Quellenhof")
9. [singularity](https://github.com/sylabs/singularity) (Tested version 3.8.5-2.el7)
10. [GATK](https://github.com/broadinstitute/gatk) (Tested version singularity image [gatk_4.4.0.0.sif](https://hub.docker.com/r/broadinstitute/gatk))
11. [DeepVariant](https://github.com/google/deepvariant) (Tested with singularity image [deepvariant_1.5.0.sif](https://hub.docker.com/r/google/deepvariant))

#### Preparing GATK singularity image


