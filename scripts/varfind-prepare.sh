#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : September 27, 2023
# Description : Filter a sample from a VCF file and create fasta file for that sample using a reference fasta
# Usage : varfind-prepare.sh -f <reference fasta,fa file> -v <VCF file> -s <sample name> -r <region>
# ./varfind-prepare.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00096
# http://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

SHORT=f:,v:,s:,r:,h
LONG=file:,name:,rename:,help
OPTS=$(getopt -a -n varfind-prepare.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: varfind-prepare.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa reference sequence file\n"
help_text+="-v | --vcf STR ground truth VCF file\n"
help_text+="-s | --sample STR Sample name to be considered from the VCF file\n"
help_text+="-r | --region STR chormosome name and region in chr:from-to format (Optional)\n"
help_text+="-h | --help Display this help message\n"

eval set -- "$OPTS"
while :
do
  case "$1" in
	-f | --file )
      file="$2"
      shift 2
      ;;
	-v | --vcf )
      vcf="$2"
      shift 2
      ;;
	-s | --sample )
      sample="$2"
      shift 2
      ;;
	-r | --region )
      region="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : varfind-prepare.sh"
      echo "Version : 1.0"
      echo -e $help_text 
      exit 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1";
	  exit 2
      ;;
  esac
done

# exit when any command fails
set -e

#Get present working directory
pwd=$(pwd)
log="${pwd}/varfinder.log"

#Function to print and log messages
function vflog(){
	echo $1;
	echo $1 >> $log
}
vflog " "
#Get Date
d=$(date)

vflog ">>> Starting varfind-prepare workflow on ${d} ..."
vflog ">>> Checking for reference fasta file ..."

if [ -z "$file" ] || [ ! -f "$file" ]; then
    vflog "The fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

file=$(realpath $file)

if [ -z "$vcf" ] || [ ! -f "$vcf" ]; then
    vflog "The vcf file ${vcf} does not exists or not specified by -v|--vcf <filename> !"
    echo -e $help_text
    exit 1;
fi

vcf=$(realpath $vcf)

if [ -z $sample ]; then
    vflog "Sample name must be specified by -s|--sample <filename> !"
    echo -e $help_text
    exit 1;
fi

vflog ">>> Checking the index file for $(basename ${file}) ..."
if [ ! -f "${file}.fai" ]; then
	vflog ">>> Index file does not exists. Creating it..."
	samtools faidx $file;
fi

vflog ">>> Checking the index file for $(basename ${vcf}) ..."
if [ ! -f "${vcf}.csi" ]; then
	vflog ">>> Index file does not exists. Creating it..."
	bcftools index $vcf;
fi

vflog ">>> Generating fasta and vcf file for the sample ${sample} ..."
if [ -z $region ]; then
	bcftools view --min-ac=1 -s $sample $vcf | bcftools norm -d all - | bgzip > "${pwd}/${sample}.vcf.gz"
	bcftools index "${pwd}/${sample}.vcf.gz"
	tabix "${pwd}/${sample}.vcf.gz"
	bcftools consensus -f $file -H A -s $sample "${pwd}/${sample}.vcf.gz" | sed "s/^>.*/>${sample}/" > "${pwd}/${sample}.fa"
	samtools faidx "${pwd}/${sample}.fa"
else
	bcftools view -r $region --min-ac=1 -s $sample $vcf | bcftools norm -d all - | bgzip > "${pwd}/${sample}.vcf.gz"
	bcftools index "${pwd}/${sample}.vcf.gz"
	tabix "${pwd}/${sample}.vcf.gz"
	samtools faidx $file $region | bcftools consensus -H A -s $sample "${pwd}/${sample}.vcf.gz" | sed "s/^>.*/>${sample}/" > "${pwd}/${sample}.fa"
	samtools faidx "${pwd}/${sample}.fa"
fi

d=$(date)
vflog ">>> Done varfind-prepare on ${d} !";