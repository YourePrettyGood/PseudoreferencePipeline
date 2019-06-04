#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <special>\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Special is a comma-separated list of options, like no_markdup or no_IR.\n"
   printf "Window size is specified bp in special, using the prefix w.\n"
   printf "Default window size is 100000.\n"
   printf "Window size of 0 indicates genome-wide average.\n"
   printf "Negative window size indicates per-scaffold averages.\n"
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
   if [[ "$i" =~ ^w[-]?[0-9][0-9]*$ ]]; then #Number *can* start with a 0
      WINDOWSIZE=${i#w} #Removes the prefixed w
   fi
done

echo "Window size selected is ${WINDOWSIZE} bp from SPECIAL string ${SPECIAL}"

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
TSVSPECIAL=""
if [[ $SPECIAL =~ "no_markdup" ]]; then
   if [[ $SPECIAL =~ "no_IR" ]]; then
      BAMSUFFIX="_sorted.bam"
      TSVSPECIAL="_nomarkdup"
   else
      BAMSUFFIX="_realigned.bam"
      TSVSPECIAL="_nomarkdup_realigned"
   fi
else
   if [[ $SPECIAL =~ "no_IR" ]]; then
      BAMSUFFIX="_sorted_markdup.bam"
   else
      BAMSUFFIX="_markdup_realigned.bam"
      TSVSPECIAL="_realigned"
   fi
fi
SAMPLEBAM="${OUTPUTDIR}${PREFIX}${BAMSUFFIX}"
FLAGSTATLOG="${OUTPUTDIR}${PREFIX}${TSVSPECIAL}_flagstat.log"

#Run samtools flagstat on the BAM:
echo "Running samtools flagstat on ${SAMPLEBAM}"
echo "${SAMTOOLS} flagstat ${SAMPLEBAM} 2>&1 > ${FLAGSTATLOG}"
${SAMTOOLS} flagstat ${SAMPLEBAM} 2>&1 > ${FLAGSTATLOG}
FSCODE=$?
if [[ $FSCODE -ne 0 ]]; then
   echo "samtools flagstat failed for ${PREFIX} on ${SAMPLEBAM} with exit code ${FSCODE}!"
   exit 2
fi
if [[ "${WINDOWSIZE}" -eq "0" ]]; then
#Specify window size of 0 for genome-wide depth:
   DEPTHTSV="${OUTPUTDIR}${PREFIX}${TSVSPECIAL}_depth_genomewide.tsv"
   SUMMARYCMD="${SCRIPTDIR}/summarizeStat.awk -v \"genomewide=1\" 2> ${OUTPUTDIR}logs/depth_summarizeStat_${PREFIX}${TSVSPECIAL}.stderr > ${DEPTHTSV}"
   SUMMARYTYPE="genome-wide"
elif [[ "${WINDOWSIZE}" -lt "0" ]]; then
#Specify negative window size for per-scaffold depth:
   DEPTHTSV="${OUTPUTDIR}${PREFIX}${TSVSPECIAL}_depth_perScaf.tsv"
   SUMMARYCMD="${SCRIPTDIR}/summarizeStat.awk 2> ${OUTPUTDIR}logs/depth_summarizeStat_${PREFIX}${TSVSPECIAL}.stderr | sort -k1,1V > ${DEPTHTSV}"
   SUMMARYTYPE="per-scaffold"
else
#Otherwise do windowed depth:
   DEPTHTSV="${OUTPUTDIR}${PREFIX}${TSVSPECIAL}_depth_w${WINDOWKB}kb.tsv"
   SUMMARYCMD="${NOW} -w ${WINDOWSIZE} -o ${DEPTHTSV} 2>&1 > ${OUTPUTDIR}logs/depth_nonOverlappingWindows_${PREFIX}${TSVSPECIAL}.log"
   SUMMARYTYPE="windowed"
fi
#Calculate windowed depth:
echo "Calculating ${SUMMARYTYPE} depth for ${SAMPLEBAM}"
DEPTHCMD="${SAMTOOLS} depth -aa ${SAMPLEBAM} 2> ${OUTPUTDIR}logs/samtoolsDepth_${PREFIX}${TSVSPECIAL}.log | ${SUMMARYCMD}"
echo "${DEPTHCMD}"
#echo "${SAMTOOLS} depth -aa ${SAMPLEBAM} 2> ${OUTPUTDIR}logs/samtoolsDepth_${PREFIX}.log | ${NOW} -w ${WINDOWSIZE} -o ${DEPTHTSV} 2>&1 > ${OUTPUTDIR}logs/depth_nonOverlappingWindows_${PREFIX}.log"
#${SAMTOOLS} depth -aa ${SAMPLEBAM} 2> ${OUTPUTDIR}logs/samtoolsDepth_${PREFIX}.log | ${NOW} -w ${WINDOWSIZE} -o ${DEPTHTSV} 2>&1 > ${OUTPUTDIR}logs/depth_nonOverlappingWindows_${PREFIX}.log
eval $DEPTHCMD
DEPTHCODE=$?
if [[ $DEPTHCODE -ne 0 ]]; then
   echo "Failed to generate ${SUMMARYTYPE} depth for ${PREFIX} on ${SAMPLEBAM} with exit code ${DEPTHCODE}!"
   exit 3
fi
echo "Flagstat and ${SUMMARYTYPE} depth finished"
