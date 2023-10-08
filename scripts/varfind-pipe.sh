#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 05, 2023
# Description : This command will execute read, map, call and compare as a pipeline
# Usage : varfind-pipe.sh -f <fasta,fa file> -n <chromosome name to filter> -r <rename chromosome (optional)>
# ./varfind-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m m -c b -t 48
# ./varfind-pipe.sh -r vgindex.giraffe.gbz -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m g -c v -t 48

SHORT=r:,f:,v:,l:,d:,s:,m:,c:,i:,t:,h
LONG=ref:,file:,vcf:,length:,depth:,sim:,mapper:,caller:,image:,threads:,help
OPTS=$(getopt -a -n varfind-pipe.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: varfind-pipe.sh [options]\n"
help_text+="-r | --ref STR reference .fasta sequence file or reference .gbz graph file\n"
help_text+="-f | --file STR .fasta file of the sample considered\n"
help_text+="-v | --vcf STR Ground truth vcf file of the sample considered\n"
help_text+="-l | --length INT read length (Default 100)\n"
help_text+="-d | --depth INT read coverage depth (Default 30)\n"
help_text+="-s | --sim STR read simulator. 'n' for NGSNGS and 'w' for wgsim (Default 'n')\n"
help_text+="-m | --mapper STR mapper/Aligner to use. 'm' for 'bwa mem', 's' for 'bwa sampe', 'b' for 'bowtie2' and 'g' for 'vg giraffe'. (Default 'm' or 'g' if the --ref is .gbz graph file)\n"
help_text+="-c | --caller STR variant caller to use. 'b' for 'bcftools', 'f' for 'freebayes', 'g' for 'gatk HaplotypeCaller', 'd' for 'DeepVariant' and 'v' for 'vg call'. (Default 'b' or 'v' if the --ref is .gbz graph file)\n"
help_text+="-i | --image STR Singularity image file if the caller run through singularity (Only for caller option 'g' and 'd')\n"
help_text+="-t | --threads INT number of threads to use (Default 'nproc')\n"
help_text+="-h | --help Display this help message\n"

eval set -- "$OPTS"
while :
do
  case "$1" in
	-r | --ref )
      ref="$2"
      shift 2
      ;;
	-f | --file )
      file="$2"
      shift 2
      ;;
	-v | --vcf )
      vcf="$2"
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
	-m | --mapper )
      mapper="$2"
      shift 2
      ;;
	-c | --caller )
      caller="$2"
      shift 2
      ;;
	-i | --image )
      image="$2"
      shift 2
      ;;
	-t | --threads )
      threads="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : varfind-pipe.sh"
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

vflog ">>> Starting varfind-pipe worlflow on ${d} ..."
vflog ">>> Checking for reference .fasta or .gbz graph file ..."

# Main file check
if [ -z "$ref" ] || [ ! -f "$ref" ]; then
    vflog "The .fasta or .gbz file ${ref} does not exists or not specified by -r|--ref <filename> !"
    echo -e $help_text
    exit 1;
fi

vflog ">>> Checking for sample sequence file ..."

if [ -z "$file" ] || [ ! -f "$file" ]; then
    vflog "The sample .fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

vflog ">>> Checking for ground truth vcf file ..."
if [ -z "$vcf" ] || [ ! -f "$vcf" ]; then
    vflog "The ground truth .vcf file ${vcf} does not exists or not specified by -v|--vcf <filename> !"
    echo -e $help_text
    exit 1;
fi

if [ ! -z "$image" ] ; then
	if [ ! -f "$image" ]; then
		vflog "The singulairty image file ${image} does not exists !"
		echo -e $help_text
		exit 1;
	fi
fi

#Read simulation parameter check
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
   vflog "No reads simulator is specified. Assigning default NGSNGS";
fi

#mapper and caller
if [ -z $threads ]; then
   threads=$(nproc);
   vflog "No number of threads are specified. Assigning default nproc=${threads}"; 
fi

if ! [[ $threads =~ $re ]] ; then
   echo "--threads must be an integer !"; exit 1;
fi

if [ ${ref##*.} == "gbz" ]; then
	vflog "The reference file is a .gbz graph"
	vflog "Assigning mapper 'vg giraffe' and ignoring any -m|--mapper option"
	mapper="g";
	vflog "Assigning caller 'vg call' and ignoring any -c|--caller option"
	caller="v";
else
	if [ -z $mapper ]; then
		vflog "No mapper is specified. Assigning default 'bwa mem'";
		mapper="m";
	fi

	if [ -z $caller ]; then
		vflog "No caller is specified. Assigning default 'bcftools call'";
		caller="b";
	fi
fi

start=`date +%s`
#Script path
sp=$(dirname "$(realpath "$0")")
prefix=${file%.*}
#Pipe START
"${sp}/varfind-reads.sh" -f $file -s $sim -l $length -d $depth
if [ ${ref##*.} == "gbz" ]; then
	"${sp}/varfind-map.sh" -g $ref -1 "${prefix}_reads_R1.fq.gz" -2 "${prefix}_reads_R2.fq.gz" -t $threads -m $mapper
	"${sp}/varfind-call.sh" -g $ref -m "${prefix}.gam" -c $caller -t $threads
else
	"${sp}/varfind-map.sh" -f $ref -1 "${prefix}_reads_R1.fq.gz" -2 "${prefix}_reads_R2.fq.gz" -t $threads -m $mapper
	if [ -z "$image" ] ; then
		"${sp}/varfind-call.sh" -f $ref -m "${prefix}.bam" -c $caller -t $threads
	else 
		"${sp}/varfind-call.sh" -f $ref -m "${prefix}.bam" -i $image -c $caller -t $threads
	fi
fi
"${sp}/varfind-compare.sh" -g $vcf -v "${prefix}.varfind.${caller}.vcf.gz" -f $file
#Pipe END
end=`date +%s`
runtime=$(($end-$start))

d=$(date)
vflog " "
vflog ">>> Parameters list"
vflog "depth,length,simulator,mapper,caller,threads,runtime"
vflog "${depth},${length},${sim},${mapper},${caller},${threads},${runtime}"
vflog ">>> Done varfind-pipe.sh on ${d} !";

