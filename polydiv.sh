#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> <special options>\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome should be line wrapped at 60 bp.\n"
   printf "Special options include markdup and IR options,\n"
   printf " variant caller, and window size (prefixed by w).\n"
   printf "Window size is in bp, used for non-overlapping windows.\n"
   exit 1
fi
#PREFIX: Prefix used for all intermediate and output files of the pipeline
PREFIX=$1
#REF: Path to the FASTA used as a reference for mapping
REF=$2
#SPECIAL: Special options, including window size, caller, and MD and IR options
SPECIAL=$3

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
fi

mkdir -p ${OUTPUTDIR}logs

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Get the MD and IR states from SPECIAL:
NOMARKDUP=""
REALIGNED=""
if [[ $SPECIAL =~ "no_markdup" ]]; then
   NOMARKDUP="_nomarkdup"
fi
if [[ ! $SPECIAL =~ "no_IR" ]]; then
   REALIGNED="_realigned"
fi
#Get the caller from SPECIAL:
CALLER=""
if [[ $SPECIAL =~ "HC" ]]; then
   CALLER="${CALLER}HC"
fi
if [[ $SPECIAL =~ "MPILEUP" ]]; then
   CALLER="${CALLER}MPILEUP"
fi

#Check that the pseudoref exists given the above options:
SAMPLE="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}_final_pseudoref.fasta"
if [[ ! -e ${SAMPLE} ]]; then
   echo "Unable to find pseudoref ${SAMPLE} given special options ${SPECIAL}"
   exit 4
fi
SAMPLEPREFIX="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}"

#Extract the window size from SPECIAL:
WINDOWSIZE=0
WINDOWSIZERE='w([0-9]+)'
if [[ ${SPECIAL} =~ $WINDOWSIZERE ]]; then
   WINDOWSIZE="${BASH_REMATCH[1]}"
   echo "Using non-overlapping windows of size ${WINDOWSIZE} bp sample ${SAMPLE}"
fi

if [[ $WINDOWSIZE -lt 1 ]]; then
   echo "Window size is less than 1, cannot proceed."
   exit 5
fi

#Note: Bash does integer arithmetic, so if you didn't use an integral
# multiple of 1000 for window size, the output filename won't be precise.
((WINDOWKB=WINDOWSIZE/1000))
POLYTSV="${SAMPLEPREFIX}_poly_w${WINDOWKB}kb.tsv"
DIVTSV="${SAMPLEPREFIX}_div_w${WINDOWKB}kb.tsv"

#Make windowed heterozygosity TSV:
${SCRIPTDIR}/oneSamplePolyDiv.sh poly ${SAMPLE} ${REF} ${WINDOWSIZE} > ${POLYTSV}
POLYCODE=$?
if [[ $POLYCODE -ne 0 ]]; then
   echo "Failed to generate windowed heterozygosity for ${SAMPLEPREFIX} with exit code ${POLYCODE}!"
   exit 2
fi
#Make windowed fixed difference TSV:
${SCRIPTDIR}/oneSamplePolyDiv.sh div ${SAMPLE} ${REF} ${WINDOWSIZE} > ${DIVTSV}
DIVCODE=$?
if [[ $DIVCODE -ne 0 ]]; then
   echo "Failed to generate windowed fixed differences for ${SAMPLE} with exit code ${DIVCODE}!"
   exit 3
fi
echo "Heterozygosity and fixed difference windows finished for ${SAMPLE}"
