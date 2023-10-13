#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 02, 2023
# Description : Construct a graph using the reference fasta file and givem samples. 
# Usage : varfind-graph.sh -f <reference fasta,fa file> -v <VCF file> -s <comma seperated sample names> -r <region>
# ./varfind-graph.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00097,HG00099,HG00100,HG00101,HG00102

SHORT=f:,v:,s:,r:,w:,h
LONG=file:,name:,rename:,write:,help
OPTS=$(getopt -a -n varfind-graph.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: varfind-prepare.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa reference sequence file\n"
help_text+="-v | --vcf STR the main VCF file\n"
help_text+="-s | --sample STR Sample names (at least 4 other than seleted sample for the simulation) to be considered from the VCF file\n"
help_text+="-r | --region STR chromosome name and region in chr:from-to format (Optional)\n"
help_text+="-w | --write STR write logs to this file (optional, default 'varfinder.log')\n"
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
	-w | --write )
      write="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : varfind-graph.sh"
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
#Log file
if [ -z "$write" ] ; then
	write='varfinder.log'
elif ! [[ $write =~ ^[0-9a-zA-Z._-]+$ ]]; then
	echo "Invalid log file name !"
	echo -e $help_text
    exit 1;
fi
log="${pwd}/${write}"

#Function to print and log messages
function vflog(){
	echo $1;
	echo $1 >> $log
}
vflog " "
#Get Date
d=$(date)

vflog ">>> Starting varfind-graph workflow on ${d} ..."
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

vflog ">>> Generating the vcf file for the samples ${sample} ..."
if [ -z $region ]; then
	bcftools view --min-ac=1 -s $sample $vcf | bcftools norm -d all - | bgzip > "${pwd}/vg.vcf.gz"
else
	bcftools view -r $region --min-ac=1 -s $sample $vcf | bcftools norm -d all - | bgzip > "${pwd}/vg.vcf.gz"
fi

vflog ">>> construct the graph and indexes ..."
tabix "${pwd}/vg.vcf.gz"
vg autoindex --workflow giraffe -r $file -v "${pwd}/vg.vcf.gz" -p vgindex

d=$(date)
vflog ">>> Done varfind-graph on ${d} !";