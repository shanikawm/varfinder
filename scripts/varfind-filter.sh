#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : September 26, 2023
# Description : This command filter a specified chromosome from a fasta file
# Usage : varfind-filter -f <fasta,fa file> -n <chromosome name to filter> -r <rename chromosome (optional)>
# ./varfind-filter.sh -f GRCh38_full_analysis_set_plus_decoy_hla.fa.gz -n chr20 -r 20

SHORT=f:,n:,r:,h
LONG=file:,name:,rename:,help
OPTS=$(getopt -a -n varfind-filter.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: varfind filter [options]\n"
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

echo ">>> Checking for fasta file ..."
if [ -z "$file" ] || [ ! -f "$file" ]; then
    >&2 echo "The fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

file=$(realpath $file)

if [ -z $name ]; then
    >&2 echo "Chromosome name must be specified by -n|--name <filename> !"
    echo -e $help_text
    exit 1;
fi

echo ">>> Checking for index file ..."
if [ ! -f "${file}.fai" ]; then
	echo ">>> Index file does not exists. Creating it..."
	samtools faidx $file;
fi

echo ">>> Filtering the Chromosome and creating the file ${name}.fa.gz ..."
if [ -z $rename ]; then
	samtools faidx $file $name | bgzip > "${pwd}/${name}.fa.gz"
else 
	samtools faidx $file $name | sed "s/^>${name}/>${rename}/" | bgzip > "${pwd}/${name}.fa.gz"
fi

echo ">>> Creating the index for ${name}.fa.gz ..."
samtools faidx "${pwd}/${name}.fa.gz"
echo "Done varfind-filter !";

