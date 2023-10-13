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
| Overall Sensitivity | 88.9661 |
| Overall Specificity | 99.9837 |
| Overall F1 Score | 87.2074 |



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


