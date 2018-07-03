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
if [[ $4 =~ "misencoded" ]]; then
   echo "Note: Attempting to convert PHRED+64 quality scores to PHRED+33"
   MISENCODED=" --fix_misencoded_quality_scores"
elif [[	$4 =~ "filter_mismatching_base_and_quals" ]]; then
   MISENCODED=" --filter_mismatching_base_and_quals"
fi
if [[ $4 =~ "no_markdup" ]]; then
   MARKDUP=""
   NOMARKDUP="_nomarkdup"
else
   MARKDUP="_markdup"
   NOMARKDUP=""
fi
if [[ $4 =~ "no_IR" ]]; then
   REALIGNED=""
else
   REALIGNED="_realigned"
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
if [[ -n "${REALIGNED}" ]]; then
   #if [[ -n "${MARKDUP}" ]]; then
      #MD IR BAM
   #else
      #noMD IR BAM
   #fi
   INPUTBAM="${OUTPUTDIR}${SAMPLE}${MARKDUP}${REALIGNED}.bam"
else
   #if [[ -n "${MARKDUP}" ]]; then
      #MD noIR BAM
   #else
      #noMD noIR BAM
   #fi
   INPUTBAM="${OUTPUTDIR}${SAMPLE}_sorted${MARKDUP}.bam"
fi
if [[ ! -e "${INPUTBAM}" ]]; then
   echo "Error: Missing input BAM ${INPUTBAM}!"
   exit 2
fi
if [[ ! -e "${INPUTBAM}.bai" ]]; then
   echo "Input BAM was missing index, creating one now."
   ${SAMTOOLS} index ${INPUTBAM}
fi

#Run HaplotypeCaller:
#Run the bamout HaplotypeCaller call if "bamout" is found in SPECIAL
echo "Starting variant calling with HaplotypeCaller for sample ${SAMPLE}"
if [[ $4 =~ "bamout" ]]; then
   java -Xmx30g -jar $GATK -T HaplotypeCaller -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct ${NPROCS} -R ${REFERENCE} -I ${INPUTBAM} -o ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_HC_ERCGVCF.vcf -bamout ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_HC_bamout.bam -forceActive -disableOptimizations${MISENCODED} 2> ${OUTPUTDIR}logs/${SAMPLE}_GATK_HaplotypeCaller.stderr > ${OUTPUTDIR}logs/${SAMPLE}_GATK_HaplotypeCaller.stdout
   HCCODE=$?
   if [[ $HCCODE -ne 0 ]]; then
      echo "GATK HaplotypeCaller on ${INPUTBAM} failed with exit code ${HCCODE}!"
      exit 3
   fi
else
   java -Xmx30g -jar $GATK -T HaplotypeCaller -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct ${NPROCS} -R ${REFERENCE} -I ${INPUTBAM} -o ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_HC_ERCGVCF.vcf${MISENCODED} 2> ${OUTPUTDIR}logs/${SAMPLE}_GATK_HaplotypeCaller.stderr > ${OUTPUTDIR}logs/${SAMPLE}_GATK_HaplotypeCaller.stdout
   HCCODE=$?
   if [[ $HCCODE -ne 0 ]]; then
      echo "GATK HaplotypeCaller on ${INPUTBAM} failed with exit code ${HCCODE}!"
      exit 3
   fi
fi
echo "HaplotypeCaller finished"

#If HaplotypeCaller successfully produced its gVCF file, run GenotypeGVCFs:
if [[ -e "${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_HC_ERCGVCF.vcf" ]]; then
   echo "Generating VCF including all sites using GenotypeGVCFs for sample ${SAMPLE}"
   java -Xmx30g -jar $GATK -T GenotypeGVCFs -nt ${NPROCS} -R ${REFERENCE} -V ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_HC_ERCGVCF.vcf -allSites -o ${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_HC_GGVCFs.vcf 2> ${OUTPUTDIR}logs/${SAMPLE}_GATK_GGVCFs.stderr > ${OUTPUTDIR}logs/${SAMPLE}_GATK_GGVCFs.stdout
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
