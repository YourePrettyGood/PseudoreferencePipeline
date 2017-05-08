#!/bin/bash
SAMPLE=$1
#SAMPLE: Prefix used for all intermediate and output files in the pipeline
REFERENCE=$2
#REFERENCE: Path to the FASTA used as a reference for mapping
SCRIPT=$3
#SCRIPT: Short name for the VCF to in.snp script to use.
#Note: It is almost universally recommended to use "Alisa"
#Note: If you append "_no_seqtk", it will skip the seqtk mutfa step.
MINDEPTH=10
#MINDEPTH: Minimum depth for filtering/masking sites (incl. invariant?)
if [[ -n $4 ]]; then
   echo "Using minimum depth of $4"
   MINDEPTH=$4
else
   echo "Using minimum depth of ${MINDEPTH}"
fi

OUTPUTDIR=""
if [[ ${SAMPLE} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${SAMPLE}`/"
   SAMPLE=`basename ${SAMPLE}`
fi

#Check that the input VCF file is there:
if [[ -e ${OUTPUTDIR}${SAMPLE}_nomarkdup_HC_GGVCFs.vcf || -e ${OUTPUTDIR}${SAMPLE}_nomarkdup_HC_GGVCFs.vcf.gz ]]; then
   echo "Auto-detected no marking of duplicates"
   NOMARKDUP="_nomarkdup"
elif [[ -e ${OUTPUTDIR}${SAMPLE}_HC_GGVCFs.vcf || -e ${OUTPUTDIR}${SAMPLE}_HC_GGVCFs.vcf.gz ]]; then
   echo "Auto-detected duplicates marked"
   NOMARKDUP=""
else
   echo "Error: Could not find input VCF."
   exit 2
fi

#Determine if input VCF is gzipped or not:
GZIPPED=""
if [[ -e ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_HC_GGVCFs.vcf.gz ]]; then
   GZIPPED=".gz"
fi
INPUTVCF="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_HC_GGVCFs.vcf${GZIPPED}"

#Load the appropriate path variables for the VCF-to-in.snp code and seqtk:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#In every case, we need to run seqtk mutfa on the in.snp to update the genome.
#Note that seqtk mutfa generates a FASTA wrapped at 60 nt per line.
if [[ ${SCRIPT} =~ "Alisa" ]]; then
   if [[ ${SCRIPT} =~ "v9" ]]; then
      echo "Running insnp_v9_alisa.py on sample ${SAMPLE} with minimum depth ${MINDEPTH}"
      ${ALISAV9} ${INPUTVCF} ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_Alisa.in.snp 20 ${MINDEPTH} 40 2 60 4 -12.5 -8.0 5
      INSNPCODE=$?
      if [[ $INSNPCODE -ne 0 ]]; then
         echo "insnp_v9_alisa.py on ${INPUTVCF} failed with exit code ${INSNPCODE}!"
         exit 3
      fi
   else
      echo "Running insnp_v8_alisa.py on sample ${SAMPLE} with minimum depth ${MINDEPTH}"
      ${ALISAV8} ${INPUTVCF} ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_Alisa.in.snp 20 ${MINDEPTH} 40 2 60 4 -12.5 -8.0 5
      INSNPCODE=$?
      if [[ $INSNPCODE -ne 0 ]]; then
         echo "insnp_v8_alisa.py on ${INPUTVCF} failed with exit code ${INSNPCODE}!"
         exit 4
      fi
   fi
   if [[ ! ${SCRIPT} =~ "_no_seqtk" ]]; then
      echo "Running seqtk mutfa on sample ${SAMPLE}"
      $SEQTK mutfa ${REFERENCE} ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_Alisa.in.snp > ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_${SCRIPT}_updated.fasta
      MUTFACODE=$?
      if [[ $MUTFACODE -ne 0 ]]; then
         echo "seqtk mutfa on ${SAMPLE}${NOMARKDUP}_Alisa.in.snp failed with exit code ${MUTFACODE}!"
         exit 5
      fi
   fi
fi
