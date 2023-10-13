# VarFinder
## Comprehensive Tool for Simulating Reads and Comparing the Read Mappings and Variant Calling Algorithms

## Installation
Make sure the dependencies are satisfied. Then, copy the files in the released script set. 
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

## Workflow
The workflow diagram with tools used and files generated in each step is shown in the below diagram. 
![varfind](https://github.com/shanikawm/varfinder/assets/8539123/d4ef6778-aef4-4f74-85d2-d8560918e0d3)

------
## Example Using Genome assembly GRCh38.p14 Chromosome 20  
### Preparing Input files

Here, we use reference [Genome assembly GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/). Download it and get the file `GCF_000001405.40_GRCh38.p14_genomic.fna`. 
Filtering the chromosome 20 can be done using the following command. 

```
./varfind-filter.sh -f GCF_000001405.40_GRCh38.p14_genomic.fna -n NC_000020.11 -r 20
```
```
./varfind-filter.sh -h
Program : varfind-filter.sh
Version : 1.0
Usage: varfind-filter.sh [options]
-f | --file STR .fasta or .fa sequence file
-n | --name STR Chromosome name to be filtered.
-r | --rename STR rename the chromosome in output sequence file (optional)
-w | --write STR write logs to this file (optional, default 'varfinder.log')
-h | --help Display this help message
```
This will produce the files `NC_000020.11.fa` and `NC_000020.11.fa.fai`
Next, we'll prepare a sequence file and a ground truth VCF file for a random sample (HG00096) of the [GRCh38](http://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz) VCF file. Download that file and use the below command. We'll restrict the example for the region 30000000-32000000

```
./varfind-prepare.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00096
```
```
./varfind-prepare.sh -h
Program : varfind-prepare.sh
Version : 1.0
Usage: varfind-prepare.sh [options]
-f | --file STR .fasta or .fa reference sequence file
-v | --vcf STR ground truth VCF file
-s | --sample STR Sample name to be considered from the VCF file
-r | --region STR chromosome name and region in chr:from-to format (Optional)
-w | --write STR write logs to this file (optional, default 'varfinder.log')
-h | --help Display this help message
```
This will produce the sequence file `HG00096.fa`, index `HG00096.fa.fai`, and the ground truth VCF file `HG00096.vcf.gz`.

If we need to consider the graph-based method as well, we have to prepare a reference graph using the below command. Here, we use five or more random samples other than the selected sample for the analysis. 

```
./varfind-graph.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00097,HG00099,HG00100,HG00101,HG00102
```
```
./varfind-graph.sh -h
Program : varfind-graph.sh
Version : 1.0
Usage: varfind-prepare.sh [options]
-f | --file STR .fasta or .fa reference sequence file
-v | --vcf STR the main VCF file
-s | --sample STR Sample names (at least 4 other than seleted sample for the simulation) to be considered from the VCF file
-r | --region STR chromosome name and region in chr:from-to format (Optional)
-w | --write STR write logs to this file (optional, default 'varfinder.log')
-h | --help Display this help message
```

### 1. Executing the workflow (ngsngs -> bwa mem -> bcftools call)
First, we'll simulate the reading using the command `varfind-reads.sh` with a read length of 100 and a coverage of 60.  

```
./varfind-reads.sh -f HG00096.fa -s n -l 100 -d 60
```
```
./varfind-reads.sh -h
Program : varfind-reads.sh
Version : 1.0
Usage: varfind-reads.sh [options]
-f | --file STR .fasta or .fa sequence file to read from
-l | --length INT read length (Default 100)
-d | --depth INT read coverage depth (Default 30)
-s | --sim STR read simulator. 'n' for NGSNGS and 'w' for wgsim (Default 'n')
-w | --write STR write logs to this file (optional, default 'varfinder.log')
-h | --help Display this help message
```

This will produce paired read files named `HG00096_reads_R1.fq.gz` and `HG00096_reads_R2.fq.gz`. 

The next step is to map the reads to the reference sequence using the command `varfind-map.sh`

```
./varfind-map.sh -f NC_000020.11.fa -1 HG00096_reads_R1.fq.gz -2 HG00096_reads_R2.fq.gz -t 48 -m m
```
```
./varfind-map.sh -h
Program : varfind-map.sh
Version : 1.0
Usage: varfind-map.sh [options]
-f | --file STR .fasta or .fa reference sequence file to map (If the mapper is not 'vg giraffe')
-g | --gbz STR reference graph .gbz file (e.g. generated by varfind-graph.sh) if the mapper is 'vg giraffe' 
-m | --mapper STR mapper/Aligner to use. 'm' for 'bwa mem', 's' for 'bwa sampe', 'b' for 'bowtie2' and 'g' for 'vg giraffe'. (Default 'm')
-1 | --read1 STR pared read .fastq or .fq file 1 
-2 | --read2 STR pared read .fastq or .fq file 2 
-t | --threads INT number of threads to use (Default 'nproc')
-w | --write STR write logs to this file (optional, default 'varfinder.log')
-h | --help Display this help message
```

This will produce a bam file called `HG00096.bam`. 

Now we can do the variant calling using `varfind-call.sh`. 

```
./varfind-call.sh -f NC_000020.11.fa -m HG00096.bam -c b -t 48
```
```
./varfind-call.sh -h
Program : varfind-call.sh
Version : 1.0
Usage: varfind-call.sh [options]
-f | --file STR .fasta or .fa reference sequence file to map (If the mapper is not 'vg call')
-g | --gbz STR reference graph .gbz file (e.g. generated by varfind-graph.sh) if the mapper is 'vg call' 
-m | --map STR .bam or .gam (if the caller is 'vg call') sequence alignment/map file
-c | --caller STR variant caller to use. 'b' for 'bcftools', 'f' for 'freebayes', 'g' for 'gatk HaplotypeCaller', 'd' for 'DeepVariant' and 'v' for 'vg call'. (Default 'b')
-i | --image STR Singularity image file if the caller run through singularity (Only for caller option 'g' and 'd')
-t | --threads INT number of threads to use (Default 'nproc')
-w | --write STR write logs to this file (optional, default 'varfinder.log')
-h | --help Display this help message
```
This step will create a VCF file called `HG00096.varfind.b.vcf.gz`. 

The final step is to compare the 2 VCF files `HG00096.varfind.b.vcf.gz` and `HG00096.vcf.gz` using `varfind-compare.sh`.

```
./varfind-compare.sh -g HG00096.vcf.gz -v HG00096.varfind.b.vcf.gz -f HG00096.fa
```
```
./varfind-compare.sh -h
Program : varfind-compare.sh
Version : 1.0
Usage: varfind-compare.sh [options]
-g | --gtvcf STR ground truth vcf file file
-v | --vfvcf STR varfinder call generated vcf file.
-f | --file STR .fasta or .fa sequence file of the sample
-w | --write STR write logs to this file (optional, default 'varfinder.log')
-h | --help Display this help message
```

This command will produce the below report with all the stats. 

| Description | Stats |
|:----------------------------------|-----------:|
| Ground Truth SNPs | 1,886 |
| Ground Truth INDELs | 271 |
| Varfind SNPs | 1,996 |
| varfind INDELs | 248 |
| SNPs Private to varfind vcf | 117 |
| INDELs Private to varfind vcf | 208 |
| Exact Matched SNPs | 1,879 |
| Exact Matched INDELs | 40 |
| True Positive (TP) | 1,919 |
| False Positive (FP) | 325 |
| True Negative (TN) | 1,997,230 |
| False Negative (FN) | 238 |
| SNP Sensitivity | 99.6288% |
| SNP Specificity | 99.9941% |
| SNP F1 Score | 96.8057% |
| INDEL Sensitivity | 14.7601% |
| INDEL Specificity | 99.9895% |
| INDEL F1 Score | 15.4142% |
| Overall Sensitivity | 88.9661% |
| Overall Specificity | 99.9837% |
| Overall F1 Score | 87.2074% |

### 2. Executing the graph-based workflow (ngsngs -> vg giraffe -> vg call)

In this scenario, the mapping step and the calling step will be different below. 

```
./varfind-map.sh -g vgindex.giraffe.gbz -1 HG00096_reads_R1.fq.gz -2 HG00096_reads_R2.fq.gz -t 48 -m g
```
This will produce a file called `HG00096.fa`, which can be used in the next step to call the variants using `varfind-call.sh`.

```
./varfind-call.sh -g vgindex.giraffe.gbz -m HG00096.gam -c v -t 48
```

### 3. Executing the whole workflow using `varfind-pipe.sh`

It is possible to execute the above-discussed workflow using the command `varfind-pipe.sh`. 
e.g:
```
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m m -c b -t 48
#For graph approach
./varfind-pipe.sh -r vgindex.giraffe.gbz -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m g -c v -t 48
```

For all 26 workflows, we can use the 26 commands. 
```
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m m -c b -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m m -c f -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m m -c g -t 48 -i ~/singularity/gatk_4.4.0.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m m -c d -t 48 -i ~/singularity/deepvariant_1.5.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m m -c b -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m m -c f -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m m -c g -t 48 -i ~/singularity/gatk_4.4.0.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m m -c d -t 48 -i ~/singularity/deepvariant_1.5.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m s -c b -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m s -c f -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m s -c g -t 48 -i ~/singularity/gatk_4.4.0.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m s -c d -t 48 -i ~/singularity/deepvariant_1.5.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m s -c b -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m s -c f -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m s -c g -t 48 -i ~/singularity/gatk_4.4.0.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m s -c d -t 48 -i ~/singularity/deepvariant_1.5.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m b -c b -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m b -c f -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m b -c g -t 48 -i ~/singularity/gatk_4.4.0.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m b -c d -t 48 -i ~/singularity/deepvariant_1.5.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m b -c b -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m b -c f -t 48
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m b -c g -t 48 -i ~/singularity/gatk_4.4.0.0.sif
./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m b -c d -t 48 -i ~/singularity/deepvariant_1.5.0.sif
./varfind-pipe.sh -r vgindex.giraffe.gbz -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m g -c v -t 48
./varfind-pipe.sh -r vgindex.giraffe.gbz -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m g -c v -t 48
```

The final report can be tabularized like below. 
| Simulator | Mapper     | Caller      | Run Time | Ground Truth SNPs | Ground Truth INDELs | Identified SNPs | Identified INDELs | Private SNPs | Private INDELs | Matched SNPs | Matched INDELs | TP   | FP    | TN      | FN   | SNP Sensitivity | SNP Specificity | SNP F1 Score | INDEL Sensitivity | INDEL Specificity | INDEL F1 Score | Overall Sensitivity | Overall Specificity | Overall F1 Score |
| :-------- | :--------- | :---------- | -------: | ----------------: | ------------------: | --------------: | ----------------: | -----------: | -------------: | -----------: | -------------: | ---: | ----: | ------: | ---: | --------------: | --------------: | -----------: | ----------------: | ----------------: | -------------: | ------------------: | ------------------: | ---------------: |
| wgsim     | bwa mem    | bcftools    | 72       | 1886              | 271                 | 1996            | 248               | 117          | 208            | 1879         | 40             | 1919 | 325   | 1997230 | 238  | 99.6288         | 99.9941         | 96.8057      | 14.7601           | 99.9895           | 15.4142        | 88.9661             | 99.9837             | 87.2074          |
| wgsim     | bwa mem    | freebayes   | 91       | 1886              | 271                 | 1876            | 239               | 144          | 239            | 1732         | 0              | 1732 | 383   | 1997172 | 425  | 91.8345         | 99.9927         | 92.0786      | 0                 | 99.988            | 0              | 80.2967             | 99.9808             | 81.0861          |
| wgsim     | bwa mem    | gatk        | 126      | 1886              | 271                 | 1885            | 254               | 21           | 21             | 1864         | 233            | 2097 | 42    | 1997513 | 60   | 98.8335         | 99.9989         | 98.8597      | 85.9778           | 99.9989           | 88.7619        | 97.2183             | 99.9978             | 97.6256          |
| wgsim     | bwa mem    | DeepVariant | 112      | 1886              | 271                 | 1943            | 259               | 73           | 27             | 1870         | 232            | 2102 | 100   | 1997455 | 55   | 99.1516         | 99.9963         | 97.6756      | 85.6088           | 99.9986           | 87.5471        | 97.4501             | 99.9949             | 96.4441          |
| ngsngs    | bwa mem    | bcftools    | 84       | 1886              | 271                 | 2072            | 7                 | 223          | 4              | 1849         | 3              | 1852 | 227   | 1997328 | 305  | 98.0381         | 99.9888         | 93.431       | 1.107             | 99.9997           | 2.1582         | 85.8599             | 99.9886             | 87.4409          |
| ngsngs    | bwa mem    | freebayes   | 107      | 1886              | 271                 | 1935            | 238               | 217          | 238            | 1718         | 0              | 1718 | 455   | 1997100 | 439  | 91.0922         | 99.9891         | 89.9241      | 0                 | 99.988            | 0              | 79.6476             | 99.9772             | 79.3533          |
| ngsngs    | bwa mem    | gatk        | 136      | 1886              | 271                 | 1772            | 230               | 12           | 17             | 1760         | 213            | 1973 | 29    | 1997526 | 184  | 93.3191         | 99.9993         | 96.2274      | 78.5977           | 99.9991           | 85.0299        | 91.4696             | 99.9985             | 94.8785          |
| ngsngs    | bwa mem    | DeepVariant | 129      | 1886              | 271                 | 2001            | 255               | 139          | 25             | 1862         | 230            | 2092 | 164   | 1997391 | 65   | 98.7274         | 99.993          | 95.8065      | 84.8708           | 99.9987           | 87.4524        | 96.9865             | 99.9917             | 94.8107          |
| wgsim     | bwa sampe  | bcftools    | 212      | 1886              | 271                 | 2053            | 242               | 173          | 202            | 1880         | 40             | 1920 | 375   | 1997180 | 237  | 99.6818         | 99.9913         | 95.4556      | 14.7601           | 99.9898           | 15.5945        | 89.0125             | 99.9812             | 86.2533          |
| wgsim     | bwa sampe  | freebayes   | 224      | 1886              | 271                 | 1968            | 232               | 238          | 232            | 1730         | 0              | 1730 | 470   | 1997085 | 427  | 91.7285         | 99.988          | 89.7768      | 0                 | 99.9883           | 0              | 80.2039             | 99.9764             | 79.4124          |
| wgsim     | bwa sampe  | gatk        | 269      | 1886              | 271                 | 1912            | 254               | 47           | 23             | 1865         | 231            | 2096 | 70    | 1997485 | 61   | 98.8865         | 99.9976         | 98.2095      | 85.2398           | 99.9988           | 88             | 97.1719             | 99.9964             | 96.9696          |
| wgsim     | bwa sampe  | DeepVariant | 249      | 1886              | 271                 | 2025            | 271               | 155          | 39             | 1870         | 232            | 2102 | 194   | 1997361 | 55   | 99.1516         | 99.9922         | 95.6277      | 85.6088           | 99.998            | 85.6088        | 97.4501             | 99.9902             | 94.4082          |
| ngsngs    | bwa sampe  | bcftools    | 213      | 1886              | 271                 | 50              | 63                | 35           | 53             | 15           | 10             | 25   | 88    | 1997467 | 2132 | 0.7953          | 99.9982         | 1.5495       | 3.69              | 99.9973           | 5.988          | 1.159               | 99.9955             | 2.2026           |
| ngsngs    | bwa sampe  | freebayes   | 240      | 1886              | 271                 | 2070            | 189               | 350          | 189            | 1720         | 0              | 1720 | 539   | 1997016 | 437  | 91.1983         | 99.9824         | 86.9565      | 0                 | 99.9905           | 0              | 79.7403             | 99.973              | 77.8985          |
| ngsngs    | bwa sampe  | gatk        | 276      | 1886              | 271                 | 1969            | 210               | 158          | 17             | 1811         | 193            | 2004 | 175   | 1997380 | 153  | 96.0233         | 99.992          | 93.9559      | 71.2177           | 99.9991           | 80.2494        | 92.9068             | 99.9912             | 92.4354          |
| ngsngs    | bwa sampe  | DeepVariant | 271      | 1886              | 271                 | 2179            | 231               | 322          | 30             | 1857         | 201            | 2058 | 352   | 1997203 | 99   | 98.4623         | 99.9838         | 91.3653      | 74.1697           | 99.9984           | 80.0796        | 95.4102             | 99.9823             | 90.1248          |
| wgsim     | bowtie2    | bcftools    | 167      | 1886              | 271                 | 4950            | 327               | 3078         | 291            | 1872         | 36             | 1908 | 3369  | 1994186 | 249  | 99.2576         | 99.8459         | 54.7688      | 13.2841           | 99.9854           | 12.0401        | 88.4561             | 99.8313             | 51.3317          |
| wgsim     | bowtie2    | freebayes   | 197      | 1886              | 271                 | 14428           | 866               | 12738        | 866            | 1690         | 0              | 1690 | 13604 | 1983951 | 467  | 89.6076         | 99.3624         | 20.7184      | 0                 | 99.9566           | 0              | 78.3495             | 99.3189             | 19.3685          |
| wgsim     | bowtie2    | gatk        | 235      | 1886              | 271                 | 3689            | 384               | 1970         | 175            | 1719         | 209            | 1928 | 2145  | 1995410 | 229  | 91.1452         | 99.9013         | 61.6681      | 77.1217           | 99.9912           | 63.8167        | 89.3834             | 99.8926             | 61.894           |
| wgsim     | bowtie2    | DeepVariant | 219      | 1886              | 271                 | 8267            | 764               | 6395         | 537            | 1872         | 227            | 2099 | 6932  | 1990623 | 58   | 99.2576         | 99.6799         | 36.8758      | 83.7638           | 99.9731           | 43.8647        | 97.311              | 99.6529             | 37.5223          |
| ngsngs    | bowtie2    | bcftools    | 182      | 1886              | 271                 | 1923            | 187               | 78           | 155            | 1845         | 32             | 1877 | 233   | 1997322 | 280  | 97.826          | 99.996          | 96.8758      | 11.8081           | 99.9922           | 13.9737        | 87.019              | 99.9883             | 87.9775          |
| ngsngs    | bowtie2    | freebayes   | 201      | 1886              | 271                 | 2381            | 226               | 656          | 226            | 1725         | 0              | 1725 | 882   | 1996673 | 432  | 91.4634         | 99.9671         | 80.853       | 0                 | 99.9886           | 0              | 79.9721             | 99.9558             | 72.4181          |
| ngsngs    | bowtie2    | gatk        | 227      | 1886              | 271                 | 1566            | 167               | 22           | 12             | 1544         | 155            | 1699 | 34    | 1997521 | 458  | 81.8663         | 99.9988         | 89.4553      | 57.1955           | 99.9993           | 70.7762        | 78.7668             | 99.9982             | 87.3521          |
| ngsngs    | bowtie2    | DeepVariant | 223      | 1886              | 271                 | 2086            | 238               | 233          | 31             | 1853         | 207            | 2060 | 264   | 1997291 | 97   | 98.2502         | 99.9883         | 93.3031      | 76.3837           | 99.9984           | 81.3359        | 95.503              | 99.9867             | 91.9437          |
| wgsim     | vg giraffe | vg call     | 254      | 1886              | 271                 | 1590            | 258               | 22           | 18             | 1568         | 240            | 1808 | 40    | 1997515 | 349  | 83.1389         | 99.9988         | 90.2186      | 88.5608           | 99.999            | 90.7372        | 83.8201             | 99.9979             | 90.2871          |
| ngsngs    | vg giraffe | vg call     | 276      | 1886              | 271                 | 1577            | 256               | 8            | 16             | 1569         | 240            | 1809 | 24    | 1997531 | 348  | 83.1919         | 99.9995         | 90.615       | 88.5608           | 99.9991           | 91.0815        | 83.8664             | 99.9987             | 90.6766          |

-------
## Appendix
#### Preparing and running GATK singularity image
```
singularity pull docker://broadinstitute/gatk:4.4.0.0
ls -sh
total 2.4G
2.4G gatk_4.4.0.0.sif
```

Running the tool

```
singularity run -B /${PWD}:/data  gatk_4.4.0.0.sif gatk

 Usage template for all tools (uses --spark-runner LOCAL when used with a Spark tool)
    gatk AnyTool toolArgs

 Usage template for Spark tools (will NOT work on non-Spark tools)
    gatk SparkTool toolArgs  [ -- --spark-runner <LOCAL | SPARK | GCS> sparkArgs ]

 Getting help
    gatk --list       Print the list of available tools

    gatk Tool --help  Print help on a particular tool

 Configuration File Specification
     --gatk-config-file                PATH/TO/GATK/PROPERTIES/FILE

 gatk forwards commands to GATK and adds some sugar for submitting spark jobs

   --spark-runner <target>    controls how spark tools are run
     valid targets are:
     LOCAL:      run using the in-memory spark runner
     SPARK:      run using spark-submit on an existing cluster 
                 --spark-master must be specified
                 --spark-submit-command may be specified to control the Spark submit command
                 arguments to spark-submit may optionally be specified after -- 
     GCS:        run using Google cloud dataproc
                 commands after the -- will be passed to dataproc
                 --cluster <your-cluster> must be specified after the --
                 spark properties and some common spark-submit parameters will be translated 
                 to dataproc equivalents

   --dry-run      may be specified to output the generated command line without running it
   --java-options 'OPTION1[ OPTION2=Y ... ]'   optional - pass the given string of options to the 
                 java JVM at runtime.  
                 Java options MUST be passed inside a single string with space-separated values.

   --debug-port <number> sets up a Java VM debug agent to listen to debugger connections on a
                         particular port number. This in turn will add the necessary java VM arguments
                         so that you don't need to explicitly indicate these using --java-options.
   --debug-suspend       sets the Java VM debug agent up so that the run get immediatelly suspended
                         waiting for a debugger to connect. By default the port number is 5005 but
                         can be customized using --debug-port

```

#### Preparing and running DeepVariant singularity image

```
singularity pull docker://google/deepvariant:1.5.0
```

Example calling variant 

```
tabix Sample_2.bam
singularity run -B ${PWD}:/data deepvariant_1.5.0.sif /opt/deepvariant/bin/run_deepvariant \
--model_type=WES --ref=/data/GRCh38_chr20.fa --reads=/data/Sample_2.bam --output_vcf=/data/dp_sample_2.vcf
```

For GPU version 

```
singularity pull docker://google/deepvariant:"1.5.0-gpu"
```


