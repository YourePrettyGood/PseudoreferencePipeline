#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <special>\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Special is a comma-separated list of options, like no_markdup or no_IR.\n"
   printf "Window size is specified bp in special, using the prefix w.\n"
   exit 1
fi
#PREFIX: Prefix used for all intermediate and output files of the pipeline
PREFIX=$1
#SPECIAL: no_markdup, no_IR, or both, plus indication of window size
SPECIAL=$2

#Default window size is 100kb:
WINDOWSIZE=100000

IFS="," read -r -a specialops <<< "${SPECIAL}"
for i in "${specialops[@]}"
   do
   if [[ "$i" =~ "^w[1-9][0-9]*$" ]]; then #Number can't start with 0
      WINDOWSIZE=${i#w} #Removes the prefixed w
   fi
done

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
if [[ $SPECIAL =~ "no_markdup" ]]; then
   if [[ $SPECIAL =~ "no_IR" ]]; then
      BAMSUFFIX="_sorted.bam"
   else
      BAMSUFFIX="_realigned.bam"
   fi
else
   if [[ $SPECIAL =~ "no_IR" ]]; then
      BAMSUFFIX="_sorted_markdup.bam"
   else
      BAMSUFFIX="_markdup_realigned.bam"
   fi
fi
SAMPLEBAM="${OUTPUTDIR}${PREFIX}${BAMSUFFIX}"
DEPTHTSV="${OUTPUTDIR}${PREFIX}_depth_w${WINDOWKB}kb.tsv"
FLAGSTATLOG="${OUTPUTDIR}${PREFIX}_flagstat.log"

#Run samtools flagstat on the BAM:
echo "Running samtools flagstat on ${SAMPLEBAM}"
echo "${SAMTOOLS} flagstat ${SAMPLEBAM} 2>&1 > ${FLAGSTATLOG}"
${SAMTOOLS} flagstat ${SAMPLEBAM} 2>&1 > ${FLAGSTATLOG}
FSCODE=$?
if [[ $FSCODE -ne 0 ]]; then
   echo "samtools flagstat failed for ${PREFIX} on ${SAMPLEBAM} with exit code ${FSCODE}!"
   exit 2
fi
#Calculate windowed depth:
echo "Calculating windowed depth for ${SAMPLEBAM}"
echo "${SAMTOOLS} depth -aa ${SAMPLEBAM} 2> ${OUTPUTDIR}logs/samtoolsDepth_${PREFIX}.log | ${NOW} -w ${WINDOWSIZE} -o ${DEPTHTSV} 2>&1 > ${OUTPUTDIR}logs/depth_nonOverlappingWindows_${PREFIX}.log"
${SAMTOOLS} depth -aa ${SAMPLEBAM} 2> ${OUTPUTDIR}logs/samtoolsDepth_${PREFIX}.log | ${NOW} -w ${WINDOWSIZE} -o ${DEPTHTSV} 2>&1 > ${OUTPUTDIR}logs/depth_nonOverlappingWindows_${PREFIX}.log
DEPTHCODE=$?
if [[ $DEPTHCODE -ne 0 ]]; then
   echo "Failed to generate windowed depth for ${PREFIX} on ${SAMPLEBAM} with exit code ${DEPTHCODE}!"
   exit 3
fi
echo "Flagstat and windowed depth finished"
