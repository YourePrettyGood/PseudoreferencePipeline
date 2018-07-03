#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <poly or div> <sample> <ref> <window size>\n"
   printf "First argument is either poly or div, to calculate heterozygosity\n"
   printf " or rate of fixed differences.\n"
   printf "Sample is a pseudoreference FASTA (optionally wrapped).\n"
   printf "Ref is the reference FASTA (identically wrapped as Sample).\n"
   printf "Window size is specified in bp, and windows are non-overlapping.\n"
   printf "Output is a TSV of scaffold, window start (1-based), and mean\n"
   printf " window value of the statistic.\n"
   printf "Sites with Ns in Sample or Ref are omitted from the calculation.\n"
   exit 1
fi
#PD: One of two strings: "poly" or "div"
# "poly" leads to calculating mean heterozygosity in windows
# "div" leads to calculating mean rate of fixed differences in windows
PD=$1
#SAMPLE: Pseudoreference FASTA wrapped at the same length as REF
# Note: SAMPLE may be unwrapped, in which case REF must also be unwrapped
SAMPLE=$2
#REF: Reference FASTA wrapped at the same length as SAMPLE
# Note: REF with IUPAC degenerate bases is not well-tested
REF=$3
#WINDOWSIZE: Non-overlapping window size in bp
WINDOWSIZE=$4

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

if [ "$PD" == "poly" ]; then
   POLYDIV="-p"
elif [ "$PD" == "div" ]; then
   POLYDIV="-d"
else
   echo "Invalid site type ${PD}"
   exit 1
fi

${LPDS} ${POLYDIV} -n ${REF} ${SAMPLE} | ${NOW} -n -w ${WINDOWSIZE}
