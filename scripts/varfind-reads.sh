#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : September 29, 2023
# Description : Silmulate reads using wgsim or ngsngs
# Usage : varfind-reads.sh -f <sample fasta,fa file> -l <read length> -d <depth> -s <simulator> 
# ./varfind-reads.sh -f HG00096.fa -s n -l 150 -d 60

SHORT=f:,l:,d:,s:,h
LONG=file:,length:,depth:,sim:,help
OPTS=$(getopt -a -n varfind-reads.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: varfind-reads.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa sequence file to read from\n"
help_text+="-l | --length INT read length (Default 100)\n"
help_text+="-d | --depth INT read coverage depth (Default 30)\n"
help_text+="-s | --sim STR read simulator. 'n' for NGSNGS and 'w' for wgsim (Default 'n')\n"
help_text+="-h | --help Display this help message\n"

eval set -- "$OPTS"
while :
do
  case "$1" in
	-f | --file )
      file="$2"
      shift 2
      ;;
	-l | --length )
      length="$2"
      shift 2
      ;;
	-d | --depth )
      depth="$2"
      shift 2
      ;;
	-s | --sim )
      sim="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : varfind-reads.sh"
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

vflog ">>> Starting varfind-reads worlflow on ${d} ..."
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

vflog ">>> Finding the length of the sample sequence file ..."
sample_len=$(awk 'BEGIN {t=0} {t+=$2} END {print t}' ${file}.fai)
vflog "Length is ${sample_len}"

re='^[0-9]+$'
if [ -z $depth ]; then
   depth=30;
   vflog "No coverage depth is specified. Assigning default 30"; 
fi
if ! [[ $depth =~ $re ]] ; then
   echo "--depth must be an integer !"; exit 1;
fi

if [ -z $length ]; then
   length=100;
   vflog "No reads length is specified. Assigning default 100";
fi

if ! [[ $length =~ $re ]] ; then
   echo "--length must be an integer !"; exit 1;
fi

if [ -z $sim ]; then
   sim="n";
fi

prefix=$(basename ${file%.*})
if [[ $sim == 'w' ]] ; then
	vflog ">>> Calculating number of reads required ..."
	reads=$(($sample_len*$depth/(2*$length))) # Divided by 2 is for paired read senario
	vflog "Need ${reads} reads for ${depth} coverage depth with ${length} reads";
	vflog ">>> Simulating reads with wgsim ..."
	wgsim -r 0 -e 0.001 -N $reads -1 $length -2 $length $file "${pwd}/${prefix}_reads_R1.fq" "${pwd}/${prefix}_reads_R2.fq";
	bgzip -f "${pwd}/${prefix}_reads_R1.fq";
	bgzip -f "${pwd}/${prefix}_reads_R2.fq";
else
	vflog ">>> Simulating reads with NGSNGS ..."
	cd $pwd;
	#ngsngs -i "../${sample}" -r $reads -l $length -seq PE -qs 20 -f fq.gz -o "${sample%.*}_reads";
	ngsngs -i $file -c $depth -l $length -seq PE -qs 30 -f fq.gz -o "${prefix}_reads";
fi

d=$(date)
vflog ">>> Done varfind-reads on ${d} !";
