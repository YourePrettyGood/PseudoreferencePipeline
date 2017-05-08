#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [number of cores] [extra option]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Number of cores allocated is 8 by default (match with --cpus-per-task).\n"
   printf "Extra option is 'misencoded' where 'misencoded' forces GATK to\n"
   printf " convert PHRED+64 qual scores to PHRED+33 qual scores.\n"
   printf "This script will automatically determine if you marked duplicates\n"
   printf " or not, and adjusts the input BAM name accordingly.\n"
   exit 1
fi
#SAMPLE: Prefix used for all intermediate and output files of the pipeline
SAMPLE=$1
#REFERENCE: Path to the FASTA used as a reference for mapping
REFERENCE=$2
#NPROCS: Number of cores to use with HaplotypeCaller
NPROCS=$3
if [[ -z "$3" ]]; then
   NPROCS=8
fi
#MISENCODED: If the BAM contains PHRED+64 quality scores, set this option
# and GATK will convert these quality scores to PHRED+33.
#Note: DO NOT set this option if you already set MISENCODED in the Indel
# Realignment step.  Only use this option if you skipped Indel Realignment.
if [[ $4 =~ "filter_misencoded_quality_scores" ]]; then
   echo "Note: Attempting to convert PHRED+64 quality scores to PHRED+33"
   MISENCODED=" --fix_misencoded_quality_scores"
elif [[	$4 =~ "filter_mismatching_base_and_quals" ]]; then
   MISENCODED=" --filter_mismatching_base_and_quals"
elif [[ -n "$4" ]]; then
   MISENCODED=" " #Non-empty value so we don't use the realigned BAM
else
   MISENCODED=""
fi

OUTPUTDIR=""
if [[ ${SAMPLE} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${SAMPLE}`/"
   SAMPLE=`basename ${SAMPLE}`
fi

mkdir -p ${OUTPUTDIR}logs

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

INPUTBAM=""
#Check for the appropriate BAM:
if [[ -n "${MISENCODED}" ]]; then
   #If we did mark duplicates, check for the markdup BAM:
   if [[ -e "${OUTPUTDIR}${SAMPLE}_sorted_markdup.bam" ]]; then
      INPUTBAM="${OUTPUTDIR}${SAMPLE}_sorted_markdup.bam"
      NOMARKDUP=""
   elif [[ -e "${OUTPUTDIR}${SAMPLE}_sorted.bam" ]]; then
      INPUTBAM="${OUTPUTDIR}${SAMPLE}_sorted.bam"
      NOMARKDUP="_nomarkdup"
      echo "Note: It appears you didn't mark duplicates for this sample."
      echo "Beware false positive variant calls."
   else
      echo "Error: Could not find appropriate input BAM file."
      exit 2
   fi
else
   #If we didn't mark duplicates for indel realignment, the post-IR BAM
   # has a slightly different name, so look for it:
   if [[ -e "${OUTPUTDIR}${SAMPLE}_realigned.bam" ]]; then
      INPUTBAM="${OUTPUTDIR}${SAMPLE}_realigned.bam"
      NOMARKDUP="_nomarkdup"
      echo "Note: It appears you didn't mark duplicates for this sample."
      echo "Beware false positive variant calls."
   elif [[ -e "${OUTPUTDIR}${SAMPLE}_markdup_realigned.bam" ]]; then
      INPUTBAM="${OUTPUTDIR}${SAMPLE}_markdup_realigned.bam"
      NOMARKDUP=""
   else
      echo "Error: Could not find appropriate input BAM file."
      exit 2
   fi
fi

#Run HaplotypeCaller:
#If a fifth option is set, output the "bamout" for diagnostics
echo "Starting variant calling with HaplotypeCaller for sample ${SAMPLE}"
if [[ -n "$5" ]]; then
   java -Xmx30g -jar $GATK -T HaplotypeCaller -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct ${NPROCS} -R ${REFERENCE} -I ${INPUTBAM} -o ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_HC_ERCGVCF.vcf -bamout ${OUTPUTDIR}${SAMPLE}_HC_bamout.bam -forceActive -disableOptimizations${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/${SAMPLE}_GATK_HaplotypeCaller.log
   HCCODE=$?
   if [[ $HCCODE -ne 0 ]]; then
      echo "GATK HaplotypeCaller on ${INPUTBAM} failed with exit code ${HCCODE}!"
      exit 3
   fi
else
   java -Xmx30g -jar $GATK -T HaplotypeCaller -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct ${NPROCS} -R ${REFERENCE} -I ${INPUTBAM} -o ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_HC_ERCGVCF.vcf${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/${SAMPLE}_GATK_HaplotypeCaller.log
   HCCODE=$?
   if [[ $HCCODE -ne 0 ]]; then
      echo "GATK HaplotypeCaller on ${INPUTBAM} failed with exit code ${HCCODE}!"
      exit 3
   fi
fi
echo "HaplotypeCaller finished"

#If HaplotypeCaller successfully produced its gVCF file, run GenotypeGVCFs:
if [[ -e "${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_HC_ERCGVCF.vcf" ]]; then
   echo "Generating VCF including all sites using GenotypeGVCFs for sample ${SAMPLE}"
   java -Xmx30g -jar $GATK -T GenotypeGVCFs -nt ${NPROCS} -R ${REFERENCE} -V ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_HC_ERCGVCF.vcf -allSites -o ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}_HC_GGVCFs.vcf 2>&1 > ${OUTPUTDIR}logs/${SAMPLE}_GATK_GGVCFs.log
   GGVCFCODE=$?
   if [[ $GGVCFCODE -ne 0 ]]; then
      echo "GATK GenotypeGVCFs on ${INPUTBAM} failed with exit code ${GGVCFCODE}!"
      exit 4
   fi
   echo "GenotypeGVCFs finished"
else
   echo "HaplotypeCaller failed for sample ${SAMPLE}"
   exit 3
fi
