#!/bin/bash
#Read in the command line arguments:
#Output usage if none supplied:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [special options]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Special options is a comma-separated list of flags to deviate\n"
   printf " from default operation.\n"
   printf "'misencoded' forces GATK to convert PHRED+64 qual scores to PHRED+33.\n"
   printf "'no_markdup' uses a BAM that doesn't have duplicates marked.\n"
   printf "'cleanup' omits primary analysis functions and instead deletes all non-log intermediate files.\n"
   exit 1
fi
#SAMPLE: Prefix used for all intermediate and output files in the pipeline
SAMPLE=$1
#REFERENCE: Path to the FASTA of the reference used for mapping
REFERENCE=$2
#SPECIAL: Special options
SPECIAL=$3
#'misencoded' forces GATK to subtract 31 from all quality scores (because
# you are declaring that all qual scores in your BAM are PHRED+64, and
# GATK requires/expects PHRED+33)
#'no_markdup' specifies that the input BAM has not had duplicates marked by
# Picard MarkDuplicates
MISENCODED=""
USEMARKDUP="_markdup"
if [[ ${SPECIAL} =~ "misencoded" ]]; then
   MISENCODED=" --fix_misencoded_quality_scores"
fi
if [[ ${SPECIAL} =~ "no_markdup" ]]; then
   USEMARKDUP=""
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

OUTPUTDIR=""
if [[ ${SAMPLE} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${SAMPLE}`/"
   SAMPLE=`basename ${SAMPLE}`
fi

mkdir -p ${OUTPUTDIR}logs

#SPECIAL options may indicate cleanup of intermediate files,
# but not log files:
if [[ $SPECIAL =~ "cleanup" ]]; then
   rm -f ${OUTPUTDIR}${SAMPLE}_realignTargets.list ${OUTPUTDIR}${SAMPLE}${USEMARKDUP}_realigned.bam ${OUTPUTDIR}${SAMPLE}${USEMARKDUP}_realigned.bai
   echo "Cleanup complete for sample ${PREFIX}"
   exit 0
fi

#Check if input files exist:
if [[ ! -e ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam ]]; then
   echo "Error: Missing input BAM file ${SAMPLE}_sorted${USEMARKDUP}.bam, cannot proceed."
   exit 2
fi
if [[ ! -e ${REFERENCE} ]]; then
   echo "Error: Missing reference FASTA file ${REFERENCE}, cannot proceed."
   exit 3
fi
if [[ ! -e ${REFERENCE}.fai ]]; then
   echo "Error: Missing reference FASTA index file ${REFERENCE}.fai, cannot proceed."
   exit 4
fi
#TODO: Add check for .dict file, since GATK is whiny about this

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Double-check the GATK .jar exists:
if [[ ! -e "${GATK}" ]]; then
   echo "GATK appears to be missing, could not find at GATK=${GATK}."
   exit 20;
fi
if [[ ! -x "$(command -v ${SAMTOOLS})" ]]; then
   echo "SAMtools appears to be missing, could not find at SAMTOOLS=${SAMTOOLS}."
   exit 18;
fi

#Make sure that the input BAM is indexed, as GATK complains about this:
if [[ ! -e ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam.bai ]]; then
   echo "Indexing input BAM ${SAMPLE}_sorted${USEMARKDUP}.bam"
   $SAMTOOLS index ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam
   BAMINDEXCODE=$?
   if [[ $BAMINDEXCODE -ne 0 ]]; then
      echo "samtools index on ${SAMPLE}_sorted${USEMARKDUP}.bam failed with exit code ${BAMINDEXCODE}!"
      exit 7
   fi
fi

#Run GATK's RealignerTargetCreator to generate a list of targets for realignment:
echo "Running GATK RealignerTargetCreator on ${SAMPLE}_sorted${USEMARKDUP}.bam"
echo "java -Xmx${JAVAMEM} -jar $GATK -T RealignerTargetCreator -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam -o ${OUTPUTDIR}${SAMPLE}_realignTargets.list${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_RTC${SAMPLE}.log"
java -Xmx${JAVAMEM} -jar $GATK -T RealignerTargetCreator -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam -o ${OUTPUTDIR}${SAMPLE}_realignTargets.list${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_RTC${SAMPLE}.log
GATKRTCCODE=$?
if [[ $GATKRTCCODE -ne 0 ]]; then
   echo "GATK RealignerTargetCreator on ${SAMPLE}_sorted${USEMARKDUP}.bam failed with exit code ${GATKRTCCODE}!"
   exit 8
fi

#If RealignerTargetCreator generated output, run IndelRealigner on that output:
if [[ -e ${OUTPUTDIR}${SAMPLE}_realignTargets.list ]]; then
   echo "Running GATK IndelRealigner on ${SAMPLE}_sorted${USEMARKDUP}.bam"
   echo "java -Xmx${JAVAMEM} -jar $GATK -T IndelRealigner -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam -targetIntervals ${OUTPUTDIR}${SAMPLE}_realignTargets.list -o ${OUTPUTDIR}${SAMPLE}${USEMARKDUP}_realigned.bam${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_IR${SAMPLE}.log"
   java -Xmx${JAVAMEM} -jar $GATK -T IndelRealigner -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam -targetIntervals ${OUTPUTDIR}${SAMPLE}_realignTargets.list -o ${OUTPUTDIR}${SAMPLE}${USEMARKDUP}_realigned.bam${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_IR${SAMPLE}.log
   GATKIRCODE=$?
   if [[ $GATKIRCODE -ne 0 ]]; then
      echo "GATK IndelRealigner on ${SAMPLE}_sorted${USEMARKDUP}.bam failed with exit code ${GATKIRCODE}!"
      exit 9
   fi
else
   echo "RealignerTargetCreator failed"
   exit 7
fi
