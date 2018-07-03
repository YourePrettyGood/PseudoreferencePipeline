#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> <window size>\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome should be line wrapped at 60 bp.\n"
   printf "Window size is in bp, used for non-overlapping windows.\n"
   exit 1
fi
#PREFIX: Prefix used for all intermediate and output files of the pipeline
PREFIX=$1
#REF: Path to the FASTA used as a reference for mapping
REF=$2
#WINDOWSIZE: Size of non-overlapping windows in bp
WINDOWSIZE=$3

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
fi

mkdir -p ${OUTPUTDIR}logs

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

((WINDOWKB=WINDOWSIZE/1000))
SAMPLE=${OUTPUTDIR}${PREFIX}*.fasta
POLYTSV=${OUTPUTDIR}${PREFIX}_poly_w${WINDOWKB}kb.tsv
DIVTSV=${OUTPUTDIR}${PREFIX}_div_w${WINDOWKB}kb.tsv

#Make windowed heterozygosity TSV:
${SCRIPTDIR}/oneSamplePolyDiv.sh poly ${SAMPLE} ${REF} ${WINDOWSIZE} > ${POLYTSV}
POLYCODE=$?
if [[ $POLYCODE -ne 0 ]]; then
   echo "Failed to generate windowed heterozygosity for ${PREFIX} with exit code ${POLYCODE}!"
   exit 2
fi
#Make windowed fixed difference TSV:
${SCRIPTDIR}/oneSamplePolyDiv.sh div ${SAMPLE} ${REF} ${WINDOWSIZE} > ${DIVTSV}
DIVCODE=$?
if [[ $DIVCODE -ne 0 ]]; then
   echo "Failed to generate windowed fixed differences for ${PREFIX} with exit code ${DIVCODE}!"
   exit 3
fi
echo "Heterozygosity and fixed difference windows finished"
