#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 04, 2023
# Description : This command will compare 2 vcf files and generate a report
# Usage : varfind-compare.sh -g <ground truth vcf file> -v <varfind generated vcf file> -f <fa file of the sample>
# ./varfind-compare.sh -g HG00096.vcf.gz -v HG00096.varfind.d.vcf.gz -f HG00096.fa


SHORT=g:,v:,f:,h
LONG=gtvcf:,vfvcf:,file:,help
OPTS=$(getopt -a -n varfind-compare.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: varfind-compare.sh [options]\n"
help_text+="-g | --gtvcf STR ground truth vcf file file\n"
help_text+="-v | --vfvcf STR varfinder call generated vcf file.\n"
help_text+="-f | --file STR .fasta or .fa sequence file of the sample\n"
help_text+="-h | --help Display this help message\n"

eval set -- "$OPTS"
while :
do
  case "$1" in
	-g | --gtvcf )
      gtvcf="$2"
      shift 2
      ;;
	-v | --vfvcf )
      vfvcf="$2"
      shift 2
      ;;
	-f | --file )
      file="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : varfind-compare.sh"
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

vflog ">>> Starting varfind-compare worlflow on ${d} ..."
vflog ">>> Checking for ground truth vcf file ..."

if [ -z "$gtvcf" ] || [ ! -f "$gtvcf" ]; then
    vflog "The ground truth vcf file ${gtvcf} does not exists or not specified by -g|--gtvcf <filename> !"
    echo -e $help_text
    exit 1;
fi

gtvcf=$(realpath $gtvcf)

vflog ">>> Checking for varfinder vcf file ..."

if [ -z "$vfvcf" ] || [ ! -f "$vfvcf" ]; then
    vflog "The varfinder vcf file ${vfvcf} does not exists or not specified by -v|--vfvcf <filename> !"
    echo -e $help_text
    exit 1;
fi

vfvcf=$(realpath $vfvcf)

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

vflog ">>> Finding ground truth SNPs and INDELs from ${vcf} ...."

snp=$(bcftools stats "${gtvcf}" | grep "number of SNPs:" | cut -f 4)
indel=$(bcftools stats "${gtvcf}" | grep "number of indels:" | cut -f 4)

vflog ">>> Comparing the ground truth vcf and varfind vcf using bcftools isec ..."
bcftools isec -c none -p varfind_compare "${gtvcf}" "${vfvcf}"

vflog ">>> Generating stats";
v_snp=$(bcftools stats "${vfvcf}" | grep "number of SNPs:" | cut -f 4)
v_indel=$(bcftools stats "${vfvcf}" | grep "number of indels:" | cut -f 4)
p_snp=$(bcftools stats varfind_compare/0001.vcf | grep "number of SNPs:" | cut -f 4)
p_indel=$(bcftools stats varfind_compare/0001.vcf | grep "number of indels:" | cut -f 4)
m_snp=$(bcftools stats varfind_compare/0002.vcf | grep "number of SNPs:" | cut -f 4)
m_indel=$(bcftools stats varfind_compare/0002.vcf | grep "number of indels:" | cut -f 4)

#INDEL Stats
indel_tp=$m_indel
indel_fp=$p_indel
indel_tn=$(($sample_len-$indel-$p_indel))
indel_fn=$(($indel-$indel_tp))
indel_sensitivity=$(bc <<< "scale=4; (${indel_tp}*100/(${indel_tp}+${indel_fn}));")
indel_specificity=$(bc <<< "scale=4; (${indel_tn}*100/(${indel_tn}+${indel_fp}));")
indel_f1=$(bc <<< "scale=4; (${indel_tp}*100/(${indel_tp}+(0.5*(${indel_fn}+${indel_fp}))));")

#SNP stats
snp_tp=$m_snp
snp_fp=$p_snp
snp_tn=$(($sample_len-$snp-$p_snp))
snp_fn=$(($snp-$snp_tp))
snp_sensitivity=$(bc <<< "scale=4; (${snp_tp}*100/(${snp_tp}+${snp_fn}));")
snp_specificity=$(bc <<< "scale=4; (${snp_tn}*100/(${snp_tn}+${snp_fp}));")
snp_f1=$(bc <<< "scale=4; (${snp_tp}*100/(${snp_tp}+(0.5*(${snp_fn}+${snp_fp}))));")

tp=$(($m_snp+$m_indel))
fp=$(($p_snp+$p_indel))
tn=$(($sample_len-$snp-$indel-$p_snp-$p_indel))
fn=$(($snp+$indel-$tp))
sensitivity=$(bc <<< "scale=4; (${tp}*100/(${tp}+${fn}));")
specificity=$(bc <<< "scale=4; (${tn}*100/(${tn}+${fp}));")
f1=$(bc <<< "scale=4; (${tp}*100/(${tp}+(0.5*(${fn}+${fp}))));")

vflog  "|  Description                      | Stats      |"
vflog  "|:----------------------------------|-----------:|"
vflog  "$(printf "|  Ground Truth SNPs                | %'10d |\n" ${snp})"                  
vflog  "$(printf "|  Ground Truth INDELs              | %'10d |\n" ${indel})"
vflog  " "                
vflog  "$(printf "|  Identified SNPs in Simulation    | %'10d |\n" ${v_snp})"              
vflog  "$(printf "|  Identified INDELs in Simulation  | %'10d |\n" ${v_indel})"
vflog  " "             
vflog  "$(printf "|  SNPs Private to Simulation       | %'10d |\n" ${p_snp})"            
vflog  "$(printf "|  INDELs Private to Simulation     | %'10d |\n" ${p_indel})"  
vflog  " "         
vflog  "$(printf "|  Exact Matched SNPs               | %'10d |\n" ${m_snp})"            
vflog  "$(printf "|  Exact Matched INDELs             | %'10d |\n" ${m_indel})"
vflog  " "             
vflog  "$(printf "|  True Positive (TP)               | %'10d |\n" ${tp})"
vflog  "$(printf "|  False Positive (FP)              | %'10d |\n" ${fp})"
vflog  "$(printf "|  True Negative (TN)               | %'10d |\n" ${tn})"
vflog  "$(printf "|  False Negative (FN)              | %'10d |\n" ${fn})"
vflog  " " 
vflog  "$(printf "|  SNP Sensitivity                  | %'9.4f%% |\n" ${snp_sensitivity})"
vflog  "$(printf "|  SNP Specificity                  | %'9.4f%% |\n" ${snp_specificity})"
vflog  "$(printf "|  SNP F1 Score                     | %'9.4f%% |\n" ${snp_f1})"
vflog  " " 
vflog  "$(printf "|  INDEL Sensitivity                | %'9.4f%% |\n" ${indel_sensitivity})"
vflog  "$(printf "|  INDEL Specificity                | %'9.4f%% |\n" ${indel_specificity})"
vflog  "$(printf "|  INDEL F1 Score                   | %'9.4f%% |\n" ${indel_f1})"
vflog  " "
vflog  "$(printf "|  ${bold}Overall Sensitivity              | %'9.4f%%${normal} |\n" ${sensitivity})" 
vflog  "$(printf "|  ${bold}Overall Specificity              | %'9.4f%%${normal} |\n" ${specificity})"
vflog  "$(printf "|  ${bold}Overall F1 Score                 | %'9.4f%%${normal} |\n" ${f1})"

d=$(date)
vflog ">>> Done varfind-compare on ${d} !";
