#!/bin/bash
#Read in the command line arguments:
#Output usage if none supplied:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <merged BAM name> [list of component BAMs]\n"
   printf "Merges and sorts BAMs using Samtools.\n"
   exit 1
fi
#MERGED: Full path and name for Merged and sorted BAM to be output
MERGED=$1
#BAMLIST: Space separated list of BAMs to merge
BAMLIST="${@:2}"

OUTPUTDIR=""
if [[ ${MERGED} =~ \/ ]]; then #If the output has a path
   OUTPUTDIR="`dirname ${MERGED}`/"
   MERGED=`basename ${MERGED}`
fi

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#if [[ ! -e ${SAMTOOLS} ]]; then
#   echo "Samtools not found, please adjust the path in the pipeline_environment script."
#   echo "pipeline_environment.sh said ${SAMTOOLS}"
#   exit 6
#fi
if [[ ! -e ${PICARD} ]]; then
   echo "Picard not found, please adjust the path in the pipeline_environment script."
   echo "pipeline_environment.sh said ${PICARD}"
   exit 6
fi

SAMPLEID=${MERGED//_sorted_markdup*bam/} #So we can merge DNA or RNA BAMs
#Adjust read groups in component BAMs to have the same SM value:
echo "Adjusting Read Groups to use SM:${SAMPLEID}"
RGBAMS=()
for file in $BAMLIST;
   do
   RGBAM=${file//.bam/_adjustedRG.bam}
   RGID=${file//_sorted_markdup*bam/}
   echo "${file} to ${RGBAM}"
#   NEWRG="@RG\\tID:${RGID}\\tLB:${RGID}\\tPL:ILLUMINA\\tSM:${SAMPLEID}"
#   ${SAMTOOLS} addreplacerg -r "${NEWRG}" ${file} > ${RGBAM}
   #Instead of samtools, use Picard to replace read groups:
   PICARDARGS="I=${file} O=${RGBAM} RGID=${RGID} RGLB=${RGID} RGPL=ILLUMINA RGPU=${RGID} RGSM=${SAMPLEID}"
   java -jar ${PICARD} AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT ${PICARDARGS}
   RGBAMS+=("${RGBAM}")
done

#Merge the input BAMs into the output BAM, maintaining sorted order:
#We assume the input BAMs are sorted by coordinate.
echo "Merging BAMs into ${MERGED}:"
echo "${BAMLIST}"
${SAMTOOLS} merge ${OUTPUTDIR}${MERGED} ${RGBAMS[@]}
MERGECODE=$?
if [[ $MERGECODE -ne 0 ]]; then
   echo "Merging of BAMs into ${MERGED} failed with exit code ${MERGECODE}!"
   exit 7
fi

#Clean up the temporary Read Group-adjusted BAMs:
rm -f ${RGBAMS}

#Index the merged BAM:
echo "Indexing merged BAM:"
${SAMTOOLS} index ${OUTPUTDIR}${MERGED}
INDEXCODE=$?
if [[ $INDEXCODE -ne 0 ]]; then
   echo "Indexing of merged BAM ${MERGED} failed with exit code ${INDEXCODE}!"
   exit 8
fi
