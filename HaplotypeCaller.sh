#!/bin/bash

#TODO: Add options for ploidy specification, probably through a
# ploidy file associated with the REF, and maybe a samples TSV
# of sample IDs and sexes.
#Is ploidy specification necessary for GATK?

#TODO: Add specification of the expected heterozygosity (-hets)
#Also for -heterozygosityStandardDeviation and -indelHeterozygosity?
#Maybe also for -stand_call_conf?

#TODO: Add specification of PCR indel model (-pcrModel) to allow for
# recommended option for PCR-free libraries

#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [number of cores] [special options]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Number of cores allocated is 8 by default (match with --cpus-per-task).\n"
   printf "Special options is a comma-separated list of flags to deviate\n"
   printf " from default operation.\n"
   printf "'misencoded' forces GATK to convert PHRED+64 qual scores to PHRED+33.\n"
   printf "'no_markdup' uses a BAM that doesn't have duplicates marked.\n"
   printf "'no_IR' uses a BAM that hasn't been indel realigned.\n"
   printf "'no_gzip' does not gzip VCFs produced by GATK (default is to gzip).\n"
   printf "'no_HC' skips the HaplotypeCaller step and proceeds to GenotypeGVCFs.\n"
   printf "'no_geno' skips the GenotypeGVCFs step.\n"
   printf "Note that the GenotypeGVCFs step does *not* perform joint genotyping.\n"
   printf "A separate task is necessary for joint genotyping.\n"
   printf "'cleanup' omits primary analysis functions and instead deletes all non-log intermediate files.\n"
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
#SPECIAL: Various special options
SPECIAL=$4
#misencoded: If the BAM contains PHRED+64 quality scores, set this option
# and GATK will convert these quality scores to PHRED+33.
#Note: DO NOT set this option if you already set MISENCODED in the Indel
# Realignment step.  Only use this option if you skipped Indel Realignment.
if [[ ${SPECIAL} =~ "misencoded" ]]; then
   echo "Note: Attempting to convert PHRED+64 quality scores to PHRED+33"
   MISENCODED=" --fix_misencoded_quality_scores"
#filter_mismatching_base_and_quals: Skips reads with differing read vs.
# qual length, though this typically only happens for truncated or
# malformatted BAMs
#In general, don't use this, as your BAM is probably wrong
elif [[	${SPECIAL} =~ "filter_mismatching_base_and_quals" ]]; then
   MISENCODED=" --filter_mismatching_base_and_quals"
fi
#no_markdup and no_IR: Specify using a BAM that hasn't been markduped or
# indel realigned
if [[ ${SPECIAL} =~ "no_markdup" ]]; then
   MARKDUP=""
   NOMARKDUP="_nomarkdup"
else
   MARKDUP="_markdup"
   NOMARKDUP=""
fi
if [[ ${SPECIAL} =~ "no_IR" ]]; then
   REALIGNED=""
else
   REALIGNED="_realigned"
fi
#no_gzip: Do not gzip the output VCFs (default is to gzip them)
if [[ ${SPECIAL} =~ "no_gzip" ]]; then
   GZSUFFIX=""
else
   GZSUFFIX=".gz"
fi
#no_HC and no_geno: Skip HaplotypeCaller or skip GenotypeGVCFs
SKIPHC=0
SKIPGENO=0
if [[ ${SPECIAL} =~ "no_HC" ]]; then
   SKIPHC=1
fi
if [[ ${SPECIAL} =~ "no_geno" ]]; then
   SKIPGENO=1
fi

#Check if memory for GATK/java is specified:
#Format of SPECIAL argument is mem_[0-9]+[MGmg]?
#The part after the underscore gets used as a suffix to -Xmx
#Since we're using a regex to capture, this is not an attack vector
JAVAMEM="30g"
JAVAMEMRE='(mem)_([0-9]+[MGmg]?)'
if [[ ${SPECIAL} =~ $JAVAMEMRE ]]; then
   JAVAMEM="${BASH_REMATCH[2]}"
fi
echo "Running java with -Xmx${JAVAMEM} for sample ${SAMPLE}"

#Extract a prefixed output directory for the sample, if provided:
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

CALLER="HC"
OUTPREFIX="${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}"
LOGPREFIX="${OUTPUTDIR}logs/${OUTPREFIX}"
INTPREFIX="${OUTPUTDIR}${OUTPREFIX}"

#SPECIAL options may indicate cleanup of intermediate files,
# but not log files:
if [[ $SPECIAL =~ "cleanup" ]]; then
   if [[ "${SKIPHC}" -eq "0" ]]; then
      rm -f ${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX} ${INTPREFIX}_bamout.bam
   fi
   if [[ "${SKIPGENO}" -eq "0" ]]; then
      rm -f ${INTPREFIX}_ERCGVCF.vcf.idx ${INTPREFIX}_GGVCFs.vcf${GZSUFFIX} ${INTPREFIX}_GGVCFs.vcf.idx
   fi
   echo "Cleanup complete for sample ${PREFIX}"
   exit 0
fi

#Check for necessary scripts/programs (GATK, SAMtools):
if [[ ! -e "${GATK}" ]]; then
   echo "GATK appears to be missing, could not find at GATK=${GATK}."
   exit 20;
fi
if [[ ! -x "$(command -v ${SAMTOOLS})" ]]; then
   echo "SAMtools appears to be missing, could not find at SAMTOOLS=${SAMTOOLS}."
   exit 18;
fi

if [[ ! -e "${INPUTBAM}" ]]; then
   echo "Error: Missing input BAM ${INPUTBAM}!"
   exit 2
fi
if [[ ! -e "${INPUTBAM}.bai" ]]; then
   echo "Input BAM was missing index, creating one now."
   ${SAMTOOLS} index ${INPUTBAM}
fi

if [[ "${SKIPHC}" -eq "0" ]]; then
   #Run HaplotypeCaller:
   #Run the bamout HaplotypeCaller call if "bamout" is found in SPECIAL
   echo "Starting variant calling with HaplotypeCaller for sample ${SAMPLE}"
   if [[ ${SPECIAL} =~ "bamout" ]]; then
      echo "java -Xmx${JAVAMEM} -jar $GATK -T HaplotypeCaller -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct ${NPROCS} -R ${REFERENCE} -I ${INPUTBAM} -o ${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX} -bamout ${INTPREFIX}_bamout.bam -forceActive -disableOptimizations${MISENCODED} 2> ${LOGPREFIX}_GATK_HaplotypeCaller.stderr > ${LOGPREFIX}_GATK_HaplotypeCaller.stdout"
      java -Xmx${JAVAMEM} -jar $GATK -T HaplotypeCaller -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct ${NPROCS} -R ${REFERENCE} -I ${INPUTBAM} -o ${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX} -bamout ${INTPREFIX}_bamout.bam -forceActive -disableOptimizations${MISENCODED} 2> ${LOGPREFIX}_GATK_HaplotypeCaller.stderr > ${LOGPREFIX}_GATK_HaplotypeCaller.stdout
      HCCODE=$?
      if [[ $HCCODE -ne 0 ]]; then
         echo "GATK HaplotypeCaller on ${INPUTBAM} failed with exit code ${HCCODE}!"
         exit 3
      fi
   else
      echo "java -Xmx${JAVAMEM} -jar $GATK -T HaplotypeCaller -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct ${NPROCS} -R ${REFERENCE} -I ${INPUTBAM} -o ${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX}${MISENCODED} 2> ${LOGPREFIX}_GATK_HaplotypeCaller.stderr > ${LOGPREFIX}_GATK_HaplotypeCaller.stdout"
      java -Xmx${JAVAMEM} -jar $GATK -T HaplotypeCaller -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct ${NPROCS} -R ${REFERENCE} -I ${INPUTBAM} -o ${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX}${MISENCODED} 2> ${LOGPREFIX}_GATK_HaplotypeCaller.stderr > ${LOGPREFIX}_GATK_HaplotypeCaller.stdout
      HCCODE=$?
      if [[ $HCCODE -ne 0 ]]; then
         echo "GATK HaplotypeCaller on ${INPUTBAM} failed with exit code ${HCCODE}!"
         exit 3
      fi
   fi
   echo "HaplotypeCaller finished for sample ${SAMPLE}"
else
   echo "Skipping HaplotypeCaller as requested by ${SPECIAL} for sample ${SAMPLE}"
fi

if [[ "${SKIPGENO}" -eq "0" ]]; then
   #If HaplotypeCaller successfully produced its gVCF file, run GenotypeGVCFs:
   if [[ -e "${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX}" ]]; then
      echo "Generating VCF including all sites using GenotypeGVCFs for sample ${SAMPLE}"
      echo "java -Xmx${JAVAMEM} -jar $GATK -T GenotypeGVCFs -nt ${NPROCS} -R ${REFERENCE} -V ${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX} -allSites -o ${INTPREFIX}_GGVCFs.vcf${GZSUFFIX} 2> ${LOGPREFIX}_GATK_GGVCFs.stderr > ${LOGPREFIX}_GATK_GGVCFs.stdout"
      java -Xmx${JAVAMEM} -jar $GATK -T GenotypeGVCFs -nt ${NPROCS} -R ${REFERENCE} -V ${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX} -allSites -o ${INTPREFIX}_GGVCFs.vcf${GZSUFFIX} 2> ${LOGPREFIX}_GATK_GGVCFs.stderr > ${LOGPREFIX}_GATK_GGVCFs.stdout
      GGVCFCODE=$?
      if [[ $GGVCFCODE -ne 0 ]]; then
         echo "GATK GenotypeGVCFs on ${INPUTBAM} failed with exit code ${GGVCFCODE}!"
         exit 4
      fi
      echo "GenotypeGVCFs finished for sample ${SAMPLE}"
   else
      echo "HaplotypeCaller failed for sample ${SAMPLE}"
      echo "${INTPREFIX}_ERCGVCF.vcf${GZSUFFIX} not found"
      exit 3
   fi
else
   echo "Skipping GenotypeGVCFs as requested by ${SPECIAL} for sample ${SAMPLE}"
fi
