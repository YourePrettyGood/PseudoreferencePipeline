#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> <special options>\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome should be line wrapped at 60 bp.\n"
   printf "Special options include markdup and IR options,\n"
   printf " variant caller, window size (prefixed by w), and\n"
   printf " jointgeno (indicating to use joint genotyping results).\n"
   printf "Window size is in bp, used for non-overlapping windows.\n"
   exit 1
fi
#PREFIX: Prefix used for all intermediate and output files of the pipeline
PREFIX=$1
#REF: Path to the FASTA used as a reference for mapping
REF=$2
#SPECIAL: Special options, including window size, caller, and MD and IR options
SPECIAL=$3

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

#Note: Bash does integer arithmetic, so if you didn't use an integral
# multiple of 1000 for window size, the output filename won't be precise.
((WINDOWKB=WINDOWSIZE/1000))
TSVSPECIAL=""
#Get the MD and IR states from SPECIAL:
NOMARKDUP=""
REALIGNED=""
if [[ $SPECIAL =~ "no_markdup" ]]; then
   NOMARKDUP="_nomarkdup"
   TSVSPECIAL="_nomarkdup"
fi
if [[ ! $SPECIAL =~ "no_IR" ]]; then
   REALIGNED="_realigned"
   TSVSPECIAL="${TSVSPECIAL}_realigned"
fi
#Get the caller from SPECIAL:
CALLER=""
if [[ $SPECIAL =~ "HC" ]]; then
   CALLER="${CALLER}HC"
fi
if [[ $SPECIAL =~ "MPILEUP" ]]; then
   if [[ ! -z ${CALLER} ]]; then
      echo "Multiple callers specified by special options ${SPECIAL}, unable to proceed."
      exit 5
   fi
   CALLER="${CALLER}MPILEUP"
fi

SAMPLEPREFIX="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}"
#vcfToPseudoref.sh simply adds _joint to the end of the prefix if
# derived from joint genotyping, so we use that here:
if [[ $SPECIAL =~ 'jointgeno' ]]; then
   SAMPLEPREFIX="${SAMPLEPREFIX}_joint"
   echo "Using jointly-genotyped results with prefix ${SAMPLEPREFIX}"
fi

#Check that the pseudoref exists given the above options:
SAMPLE="${SAMPLEPREFIX}_final_pseudoref.fasta"
if [[ ! -e ${SAMPLE} ]]; then
   echo "Unable to find pseudoref ${SAMPLE} given special options ${SPECIAL}"
   exit 4
fi

#Check that necessary scripts/programs exist:
if [[ ! -x "$(command -v ${NOW})" ]]; then
   echo "nonOverlappingWindows appears to be missing, could not find at NOW=${NOW}."
   exit 22;
fi
if [[ ! -e "${SCRIPTDIR}/summarizeStat.awk" ]]; then
   echo "summarizeStat.awk appears to be missing from the installation at ${SCRIPTDIR}/summarizeStat.awk.";
   exit 23;
fi
if [[ ! -x "$(command -v ${LPDS})" ]]; then
   echo "listPolyDivSites appears to be missing, could not find at LPDS=${LPDS}."
   exit 24;
fi

if [[ "${WINDOWSIZE}" -eq "0" ]]; then
#Specify window size of 0 for genome-wide poly and div:
   POLYTSV="${SAMPLEPREFIX}_poly_genomewide.tsv"
   DIVTSV="${SAMPLEPREFIX}_div_genomewide.tsv"
   SUMMARYCMD="${SCRIPTDIR}/summarizeStat.awk -v \"genomewide=1\""
   SUMMARYTYPE="genome-wide"
elif [[ "${WINDOWSIZE}" -lt "0" ]]; then
#Specify negative window size for per-scaffold poly and div:
   POLYTSV="${SAMPLEPREFIX}_poly_w${WINDOWKB}kb.tsv"
   DIVTSV="${SAMPLEPREFIX}_div_w${WINDOWKB}kb.tsv"
   SUMMARYCMD="${SCRIPTDIR}/summarizeStat.awk"
   SUMMARYTYPE="per-scaffold"
elif [[ "${WINDOWSIZE}" -eq "1" ]]; then
   POLYTSV="${SAMPLEPREFIX}_poly_persite.tsv.gz"
   DIVTSV="${SAMPLEPREFIX}_div_persite.tsv.gz"
   SUMMARYCMD="gzip -9"
   SUMMARYTYPE="per-site"
else
#Otherwise do windowed poly and div:
   POLYTSV="${SAMPLEPREFIX}_poly_w${WINDOWKB}kb.tsv"
   DIVTSV="${SAMPLEPREFIX}_div_w${WINDOWKB}kb.tsv"
   SUMMARYCMD="${NOW} -n -w ${WINDOWSIZE}"
   SUMMARYTYPE="windowed"
fi

#Make windowed heterozygosity TSV:
echo "Calculating ${SUMMARYTYPE} heterozygosity for ${SAMPLE}"
POLYCMD="${LPDS} -p -n ${REF} ${SAMPLE} | ${SUMMARYCMD} > ${POLYTSV}"
eval $POLYCMD
#${SCRIPTDIR}/oneSamplePolyDiv.sh poly ${SAMPLE} ${REF} ${WINDOWSIZE} > ${POLYTSV}
POLYCODE=$?
if [[ $POLYCODE -ne 0 ]]; then
   echo "Failed to generate heterozygosity for ${SAMPLEPREFIX} with exit code ${POLYCODE}!"
   exit 2
fi
#Make windowed fixed difference TSV:
echo "Calculating ${SUMMARYTYPE} fixed difference rate for ${SAMPLE}"
DIVCMD="${LPDS} -d -n ${REF} ${SAMPLE} | ${SUMMARYCMD} > ${DIVTSV}"
eval $DIVCMD
#${SCRIPTDIR}/oneSamplePolyDiv.sh div ${SAMPLE} ${REF} ${WINDOWSIZE} > ${DIVTSV}
DIVCODE=$?
if [[ $DIVCODE -ne 0 ]]; then
   echo "Failed to generate fixed differences for ${SAMPLE} with exit code ${DIVCODE}!"
   exit 3
fi
echo "Heterozygosity and fixed difference windows finished for ${SAMPLE}"
