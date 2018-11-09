#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [number of cores] [extra option]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Number of cores allocated is 8 by default (match with --cpus-per-task).\n"
   printf "Extra option is for indicating if you marked duplicates and/or\n"
   printf "realigned with GATK IndelRealigner\n"
   printf "Signify using the BAM without duplicates marked with no_markdup\n"
   printf "Signify using a non-indel-realigned BAM with no_IR\n"
   printf "Separate multiple extra options with commas and no spaces\n"
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

if [[ $4 =~ "no_markdup" ]]; then
   MARKDUP=""
   NOMARKDUP="_nomarkdup"
else
   MARKDUP="_markdup"
   NOMARKDUP=""
fi
if [[ $4 =~ "no_IR" ]]; then
   REALIGNED=""
else
   REALIGNED="_realigned"
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

BCFTOOLSVERSION="`${BCFTOOLS} --version | awk '/^bcftools/{split($2, versionarr, "-"); print versionarr[1];}'`"
VERSIONARR=(${BCFTOOLSVERSION//./ })
$SAMTOOLS mpileup -ugf ${REF} ${INPUTBAM} 2> ${OUTPUTDIR}logs/samtoolsMpileup${PREFIX}.stderr
| $BCFTOOLS call -m -Oz -o ${OUTPUTVCF} 2>&1 > ${OUTPUTDIR}logs/bcftoolsCall${PREFIX}.log
if [[ "${VERSIONARR[0]}" -gt "0" ]]; then
   if [[ "${VERSIONARR[0]}" -eq "1" && "${VERSIONARR[1]}" -lt "7" ]]; then
      echo "Using samtools mpileup, since BCFtools version is < 1.7"
      MPILEUPCMD="${SAMTOOLS} mpileup -ugf ${REF}"
      CALLCMD="${BCFTOOLS} call -m -Oz -o"
   else
      echo "Using bcftools mpileup, since BCFtools version is >= 1.7"
      MPILEUPCMD="${BCFTOOLS} mpileup -Ou -o - --threads ${NUMPROCS} -f ${REF}"
      CALLCMD="${BCFTOOLS} call --threads ${NUMPROCS} -m -Oz -o"
   fi
else
   echo "Your version of BCFtools (${BCFTOOLSVERSION}) is too old to use, please update it"
   exit 8
fi

INPUTBAM=""
#Check for the appropriate BAM:
if [[ -n "${REALIGNED}" ]]; then
   #if [[ -n "${MARKDUP}" ]]; then
      #MD IR BAM
   #else
      #noMD IR BAM
   #fi
   INPUTBAM="${OUTPUTDIR}${PREFIX}${MARKDUP}${REALIGNED}.bam"
else
   #if [[ -n "${MARKDUP}" ]]; then
      #MD noIR BAM
   #else
      #noMD noIR BAM
   #fi
   INPUTBAM="${OUTPUTDIR}${PREFIX}_sorted${MARKDUP}.bam"
fi
if [[ ! -e "${INPUTBAM}" ]]; then
   echo "Error: Missing input BAM ${INPUTBAM}!"
   exit 2
fi

OUTPUTVCF="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_mpileupcall.vcf.gz"
#Detect variants (including indels) with mpileup, and call variants with bcftools call, storing as a bgzipped VCF:
${MPILEUPCMD} ${INPUTBAM} 2> ${OUTPUTDIR}logs/samtoolsMpileup${PREFIX}.stderr | ${CALLCMD} ${OUTPUTVCF} 2>&1 > ${OUTPUTDIR}logs/bcftoolsCall${PREFIX}.log
MPCALLCODE=$?
if [[ $MPCALLCODE -ne 0 ]]; then
   echo "Samtools Mpileup or BCFtools call failed on ${INPUTBAM} with exit code ${MPCALLCODE}!"
   exit 3
fi
#Index the bgzipped VCF:
${TABIX} ${OUTPUTVCF}
TABIXCODE=$?
if [[ $TABIXCODE -ne 0 ]]; then
   echo "Tabix failed to index ${OUTPUTVCF} with exit code ${TABIXCODE}!"
   exit 4
fi
echo "Samtools variant calling finished"
