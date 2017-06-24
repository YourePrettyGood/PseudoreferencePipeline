#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [number of cores] [extra option]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Number of cores allocated is 8 by default (match with --cpus-per-task).\n"
   printf "Extra option is just if you realigned with GATK IndelRealigner\n"
   printf "This script will automatically determine if you marked duplicates\n"
   printf " or not, and adjusts the input BAM name accordingly.\n"
   exit 1
fi
#PREFIX: Prefix used for all intermediate and output files of the pipeline
PREFIX=$1
#REF: Path to the FASTA used as a reference for mapping
REF=$2
#NUMPROCS: Number of cores to use with samtools and bcftools
NUMPROCS=$3
if [[ -z "$3" ]]; then
   NUMPROCS=8
fi

#REALIGNED: Whether to use the indel-realigned BAM for variant calling
if [[ -n "$4" ]]; then
   REALIGNED="realigned" #Non-empty value so we do use the realigned BAM
else
   REALIGNED="" #Empty value to avoid the realigned BAM
fi

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
fi

mkdir -p ${OUTPUTDIR}logs

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

INPUTBAM=""
#Check for the appropriate BAM:
if [[ -n "${REALIGNED}" ]]; then
   #If we didn't mark duplicates for indel realignment, the post-IR BAM
   # has a slightly different name, so look for it:
   if [[ -e "${OUTPUTDIR}${PREFIX}_realigned.bam" ]]; then
      INPUTBAM="${OUTPUTDIR}${PREFIX}_realigned.bam"
      NOMARKDUP="_nomarkdup"
      echo "Note: It appears you didn't mark duplicates for this sample."
      echo "Beware false positive variant calls."
   elif [[ -e "${OUTPUTDIR}${PREFIX}_markdup_realigned.bam" ]]; then
      INPUTBAM="${OUTPUTDIR}${PREFIX}_markdup_realigned.bam"
      NOMARKDUP=""
   else
      echo "Error: Could not find appropriate input BAM file."
      exit 2
   fi
else
   #If we did mark duplicates, check for the markdup BAM:
   if [[ -e "${OUTPUTDIR}${PREFIX}_sorted_markdup.bam" ]]; then
      INPUTBAM="${OUTPUTDIR}${PREFIX}_sorted_markdup.bam"
      NOMARKDUP=""
   elif [[ -e "${OUTPUTDIR}${PREFIX}_sorted.bam" ]]; then
      INPUTBAM="${OUTPUTDIR}${PREFIX}_sorted.bam"
      NOMARKDUP="_nomarkdup"
      echo "Note: It appears you didn't mark duplicates for this sample."
      echo "Beware false positive variant calls."
   else
      echo "Error: Could not find appropriate input BAM file."
      exit 2
   fi
fi

OUTPUTVCF="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_mpileupcall.vcf.gz"
#Detect variants (including indels) with mpileup, and call variants with bcftools call, storing as a bgzipped VCF:
$SAMTOOLS mpileup -ugf ${REF} ${INPUTBAM} 2> ${OUTPUTDIR}logs/samtoolsMpileup${PREFIX}.stderr | $BCFTOOLS call -m -Oz -o ${OUTPUTVCF} 2>&1 > ${OUTPUTDIR}logs/bcftoolsCall${PREFIX}.log
MPCALLCODE=$?
if [[ $MPCALLCODE -ne 0 ]]; then
   echo "Samtools Mpileup or BCFtools call failed on ${INPUTBAM} with exit code ${MPCALLCODE}!"
   exit 3
fi
#Index the bgzipped VCF:
${TABIX} ${OUTPUTVCF}
TABIXCODE=$?
if [[ $TABIXCODE -ne 0 ]]; then
   echo "Tabix failed to index ${INPUTBAM} with exit code ${TABIXCODE}!"
   exit 4
fi
echo "Samtools variant calling finished"
