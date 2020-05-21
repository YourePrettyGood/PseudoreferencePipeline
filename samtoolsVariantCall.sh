#!/bin/bash

#TODO: Add options for ploidy specification, probably through a
# ploidy file associated with the REF, and maybe a samples TSV
# of sample IDs and sexes.
#Relevant questions:
#1) How does the user specify the ploidy file?
#2) How does the user specify the samples TSV?
#Potential answers:
#1.1) ${REF}.ploidy file, and we auto-detect existence to add it to the call
#2.1) ${METADATA}.sexes file, and we auto-detect existence to add it to the call
#2.1) The problem here is that ${METADATA} isn't passed, so we lack that info
#2.2) Append to ${SPECIAL}, preferably automatically from *ArrayCall_v2.sh
#2.2) Could make the option something like 'sex_[MF]', easy to parse
#For now, we simply use the defaults where everything is assumed diploid.

#TODO: Add specification of the "expected substitution rate" (-P)

#TODO: Add -f GQ for additional GQ annotation to bcftools call

#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [number of cores] [special options]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Number of cores allocated is 8 by default (match with --cpus-per-task).\n"
   printf "Special options is a comma-separated list of flags to deviate\n"
   printf " from default operation.\n"
   printf "'no_markdup' indicates using a BAM without duplicates marked as input.\n"
   printf "'no_IR' indicates using a BAM without realigning around indels as\n"
   printf " input.\n"
   printf "'hets_#' gives the heterozygosity/substitution rate to use.\n"
   printf " This is generally your expected theta, default is to skip this.\n"
   printf "'cleanup' omits primary analysis functions and instead deletes all non-log intermediate files.\n"
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
#SPECIAL: Options/flags for non-default operation
SPECIAL=$4
if [[ ${SPECIAL} =~ "no_markdup" ]]; then
   MARKDUP=""
   NOMARKDUP="_nomarkdup"
else
   MARKDUP="_markdup"
   NOMARKDUP=""
fi
if [[ ${SPECIAL} =~ "no_IR" ]]; then
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

OUTPUTVCF="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_mpileupcall.vcf.gz"

#SPECIAL options may indicate cleanup of intermediate files,
# but not log files:
if [[ $SPECIAL =~ "cleanup" ]]; then
   rm -f ${OUTPUTVCF} ${OUTPUTVCF}.tbi
   echo "Cleanup complete for sample ${PREFIX}"
   exit 0
fi

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Check if expected heterozygosity is specified:
#Format of SPECIAL argument is hets_[01]?[.]?[0-9]*
THETA=""
THETARE='(hets)_([01]?[.]?[0-9]*)'
if [[ ${SPECIAL} =~ $THETARE ]]; then
   THETA="${BASH_REMATCH[2]}"
fi

#Check that the necessary scripts/tools exist (excluding awk and java):
if [[ ! -x "$(command -v ${BCFTOOLS})" ]]; then
   echo "BCFtools appears to be missing, could not find at BCFTOOLS=${BCFTOOLS}."
   exit 17;
fi
if [[ ! -x "$(command -v ${SAMTOOLS})" ]]; then
   echo "SAMtools appears to be missing, could not find at SAMTOOLS=${SAMTOOLS}."
   exit 18;
fi
if [[ ! -x "$(command -v ${TABIX})" ]]; then
   echo "tabix appears to be missing, could not find at TABIX=${TABIX}."
   exit 19;
fi

#Add in the expected heterozygosity argument if specified:
THETAOPT=""
if [[ -n "${THETA}" ]]; then
   THETAOPT="-P ${THETA}"
fi

BCFTOOLSVERSION=`${BCFTOOLS} --version | awk '/^bcftools/{split($2, versionarr, "-"); print versionarr[1];}'`
VERSIONARR=(${BCFTOOLSVERSION//./ })
if [[ "${VERSIONARR[0]}" -gt "0" ]]; then
   if [[ "${VERSIONARR[0]}" -eq "1" && "${VERSIONARR[1]}" -lt "7" ]]; then
      echo "Using samtools mpileup, since BCFtools version is < 1.7"
      if [[ "${VERSIONARR[1]}" -lt "3" ]]; then
         ANNOTATIONS="DP,DP4,SP"
      else
         ANNOTATIONS="DP,AD,SP,ADF,ADR"
      fi
      MPILEUPCMD="${SAMTOOLS} mpileup -ugf ${REF} -t ${ANNOTATIONS}"
      CALLCMD="${BCFTOOLS} call -m ${THETAOPT} -Oz -o"
   else
      echo "Using bcftools mpileup, since BCFtools version is >= 1.7"
      ANNOTATIONS="FORMAT/DP,FORMAT/AD,FORMAT/SP,FORMAT/ADF,FORMAT/ADR"
      MPILEUPCMD="${BCFTOOLS} mpileup -Ou -o - --threads ${NUMPROCS} -f ${REF} -a ${ANNOTATIONS}"
      CALLCMD="${BCFTOOLS} call --threads ${NUMPROCS} -m ${THETAOPT} -Oz -o"
   fi
else
   echo "Your version of BCFtools (${BCFTOOLSVERSION}) is too old to use, please update it"
   exit 8
fi

INPUTBAM=""
#Check for the appropriate BAM:
if [[ -n "${REALIGNED}" ]]; then
   INPUTBAM="${OUTPUTDIR}${PREFIX}${MARKDUP}${REALIGNED}.bam"
else
   INPUTBAM="${OUTPUTDIR}${PREFIX}_sorted${MARKDUP}.bam"
fi
if [[ ! -e "${INPUTBAM}" ]]; then
   echo "Error: Missing input BAM ${INPUTBAM}!"
   exit 2
fi

#Detect variants (including indels) with mpileup, and call variants with bcftools call, storing as a bgzipped VCF:
echo "Calling variants with MPILEUP for sample ${PREFIX}"
echo "${MPILEUPCMD} ${INPUTBAM} 2> ${OUTPUTDIR}logs/samtoolsMpileup${PREFIX}.stderr | ${CALLCMD} ${OUTPUTVCF} 2>&1 > ${OUTPUTDIR}logs/bcftoolsCall${PREFIX}.log"
${MPILEUPCMD} ${INPUTBAM} 2> ${OUTPUTDIR}logs/samtoolsMpileup${PREFIX}.stderr | ${CALLCMD} ${OUTPUTVCF} 2>&1 > ${OUTPUTDIR}logs/bcftoolsCall${PREFIX}.log
MPCALLCODE=$?
if [[ $MPCALLCODE -ne 0 ]]; then
   echo "Samtools Mpileup or BCFtools call failed on ${INPUTBAM} with exit code ${MPCALLCODE}!"
   exit 3
fi
#Index the bgzipped VCF:
echo "Indexing the MPILEUP VCF for sample ${PREFIX}"
echo "${TABIX} ${OUTPUTVCF}"
${TABIX} ${OUTPUTVCF}
TABIXCODE=$?
if [[ $TABIXCODE -ne 0 ]]; then
   echo "Tabix failed to index ${OUTPUTVCF} with exit code ${TABIXCODE}!"
   exit 4
fi
echo "Samtools variant calling finished for sample ${PREFIX}"
