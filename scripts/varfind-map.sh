#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : September 29, 2023
# Description : Align reads to a reference
# Usage : varfind-map.sh -f <reference fasta,fa file> -m <mapper> -r <read fq file 1> -r <read fq file 2> -t <no of threads>
# /varfind-map.sh -f chr20.fa.gz -1 HG00096.fasta_reads_R1.fq.gz -2 HG00096.fasta_reads_R2.fq.gz -t 48

SHORT=f:,m:,1:,2:,t:,h
LONG=file:,mapper:,read1:,read2:,help
OPTS=$(getopt -a -n varfind-map.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: varfind-map.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa reference sequence file to map\n"
help_text+="-m | --mapper STR mapper/Aligner to use. 'm' for 'bwa mem', 's' for 'bwa sampe', 'b' for 'bowtie2' and 'g' for 'vg giraffe'. (Default 'm')\n"
help_text+="-1 | --read1 STR pared read .fastq or .fq file 1 \n"
help_text+="-2 | --read2 STR pared read .fastq or .fq file 2 \n"
help_text+="-t | --threads INT number of threads to use (Default 'nproc')\n"
help_text+="-h | --help Display this help message\n"

eval set -- "$OPTS"
while :
do
  case "$1" in
	-f | --file )
      file="$2"
      shift 2
      ;;
	-m | --mapper )
      mapper="$2"
      shift 2
      ;;
	-1 | --read1 )
      read1="$2"
      shift 2
      ;;
	-2 | --read2 )
      read2="$2"
      shift 2
      ;;
	-t | --threads )
      threads="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : varfind-map.sh"
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

vflog ">>> Starting varfind-map worlflow on ${d} ..."
vflog ">>> Checking for sequence fasta file ..."

if [ -z "$file" ] || [ ! -f "$file" ]; then
    vflog "The fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

file=$(realpath $file)

vflog ">>> Checking the index file for $(basename ${file}) ..."
if [ ! -f "${file}.fai" ]; then
	vflog ">>> Index file does not exists. Creating it..."
	samtools faidx $file;
fi

if [ -z "$read1" ] || [ ! -f "$read1" ]; then
    vflog "The fastq file 1 ${read1} does not exists or not specified by -1|--read1 <filename> !"
    echo -e $help_text
    exit 1;
fi

read1=$(realpath $read1)

if [ -z "$read2" ] || [ ! -f "$read2" ]; then
    vflog "The fastq file 2 ${read2} does not exists or not specified by -2|--read2 <filename> !"
    echo -e $help_text
    exit 1;
fi

read2=$(realpath $read2)

if [ -z $mapper ]; then
   mapper="m";
fi

re='^[0-9]+$'
if [ -z $threads ]; then
   threads=$(nproc);
   vflog "No number of threads are specified. Assigning default nproc=${threads}"; 
fi
if ! [[ $threads =~ $re ]] ; then
   echo "--threads must be an integer !"; exit 1;
fi

prefix=$(basename ${read1%.*.*.*})
if [[ $mapper == 'm' ]] ; then
	vflog ">>> Aligning reads to the reference \"$(basename $file)\" using \"bwa mem\"";
	bwa mem -t $threads -R "@RG\tID:${prefix}\tSM:${prefix}\tLB:L1" "${file}" "${read1}" "${read2}" | samtools sort - > "${pwd}/${prefix}.bam"
elif [[ $mapper == 's' ]] ; then
	vflog ">>> Aligning reads to the reference \"$(basename $file)\" using \"bwa sampe\"";
	bwa aln "${file}" "${read1}" > "${read1}.sai"
	bwa aln "${file}" "${read2}" > "${read2}.sai"
	bwa sampe -r "@RG\tID:${prefix}\tSM:${prefix}\tLB:L1" "${file}" "${read1}.sai" "${read2}.sai" "${read1}" "${read2}" | samtools sort - > "${pwd}/${prefix}.bam"
elif [[ $mapper == 'b' ]] ; then
	vflog ">>> Aligning reads to the reference \"$(basename $file)\" using \"bowtie2\"";
	bowtie2-build "${file}" "${file}.bt2-idx"
	bowtie2 -x "${file}.bt2-idx" -p $threads -1 "${read1}" -2 "${read2}" | samtools sort - > "${pwd}/${prefix}.bam"
else
	vflog "Unknown option \"${mapper}\" for --mapper!"
    echo -e $help_text
    exit 1;
fi

d=$(date)
vflog ">>> Done varfind-map on ${d} !";
