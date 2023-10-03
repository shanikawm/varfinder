#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : September 26, 2023
# Description : This command filter a specified chromosome from a fasta file
# Usage : varfind-filter -f <fasta,fa file> -n <chromosome name to filter> -r <rename chromosome (optional)>
# ./varfind-filter.sh -f GRCh38_full_analysis_set_plus_decoy_hla.fa.gz -n chr20 -r 20
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
# GCF_000001405.40_GRCh38.p14_genomic.fna -> NC_000020.11 Homo sapiens chromosome 20, GRCh38.p14 Primary Assembly
# ./varfind-filter.sh -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz -n NC_000020.11 -r 20

SHORT=f:,n:,r:,h
LONG=file:,name:,rename:,help
OPTS=$(getopt -a -n varfind-filter.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: varfind-filter.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa sequence file\n"
help_text+="-n | --name STR Chromosome name to be filtered.\n"
help_text+="-r | --rename STR rename the chromosome in output sequence file (optional)\n"
help_text+="-h | --help Display this help message\n"

eval set -- "$OPTS"
while :
do
  case "$1" in
	-f | --file )
      file="$2"
      shift 2
      ;;
	-n | --name )
      name="$2"
      shift 2
      ;;
	-r | --rename )
      rename="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : varfind-filter.sh"
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

vflog ">>> Starting varfind-filter worlflow on ${d} ..."
vflog ">>> Checking for fasta file ..."

if [ -z "$file" ] || [ ! -f "$file" ]; then
    vflog "The fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

file=$(realpath $file)

if [ -z $name ]; then
    vflog "Chromosome name must be specified by -n|--name <filename> !"
    echo -e $help_text
    exit 1;
fi

vflog ">>> Checking for index file ..."
if [ ! -f "${file}.fai" ]; then
	vflog ">>> Index file does not exists. Creating it..."
	samtools faidx $file;
fi

vflog ">>> Filtering the Chromosome and creating the file ${name}.fa ..."
if [ -z $rename ]; then
	samtools faidx $file $name > "${pwd}/${name}.fa"
else 
	samtools faidx $file $name | sed "s/^>${name}/>${rename}/" > "${pwd}/${name}.fa"
fi

vflog ">>> Creating the index for ${name}.fa ..."
samtools faidx "${pwd}/${name}.fa"
bwa index "${pwd}/${name}.fa";
d=$(date)
vflog ">>> Done varfind-filter on ${d} !";

