#!/bin/bash
#Read in the command line arguments:
#Output usage if none supplied:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [extra option]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Extra option is either 'misencoded' or 'no_markdup'\n"
   printf "where 'misencoded' forces GATK to convert PHRED+64 qual scores to\n"
   printf "PHRED+33 qual scores, and 'no_markdup' uses the pre-markdup\n"
   printf "BAM for indel realignment.\n"
   exit 1
fi
#SAMPLE: Prefix used for all intermediate and output files in the pipeline
SAMPLE=$1
#REFERENCE: Path to the FASTA of the reference used for mapping
REFERENCE=$2
#MISENCODED: Flag used if your BAM has PHRED+64 quality scores, as GATK will
# complain about this.  If you put anything in MISENCODED, it will trigger
# usage of --fix_misencoded_quality_scores, which makes GATK subtract 31 from
# all quality scores.
#MISENCODED: Flag used to trigger conversion of PHRED+64 quality scores in
# your BAM to PHRED+33, because GATK does not auto-convert PHRED+64, and
# complains if you give it PHRED+64 data...
#USEMARKDUP: Flag used to trigger use of non-markdupped BAM for Indel
# Realignment
if [[ -z "$3" ]]; then
   MISENCODED=""
   USEMARKDUP="_markdup"
elif [[ $3 =~ "no_markdup" ]]; then
   USEMARKDUP=""
   MISENCODED=""
else
   USEMARKDUP="_markdup"
   MISENCODED=" --fix_misencoded_quality_scores"
fi

OUTPUTDIR=""
if [[ ${SAMPLE} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${SAMPLE}`/"
   SAMPLE=`basename ${SAMPLE}`
fi

mkdir -p ${OUTPUTDIR}logs

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
if [[ ! -e ${GATK} ]]; then
   echo "GATK not found, please adjust the path in the pipeline_environment script."
   echo "pipeline_environment.sh said ${GATK}"
   exit 5
fi
if [[ ! -e ${SAMTOOLS} ]]; then
   echo "Samtools not found, please adjust the path in the pipeline_environment script."
   echo "pipeline_environment.sh said ${SAMTOOLS}"
   exit 6
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
java -Xmx30g -jar $GATK -T RealignerTargetCreator -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam -o ${OUTPUTDIR}${SAMPLE}_realignTargets.list${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_RTC${SAMPLE}.log
GATKRTCCODE=$?
if [[ $GATKRTCCODE -ne 0 ]]; then
   echo "GATK RealignerTargetCreator on ${SAMPLE}_sorted${USEMARKDUP}.bam failed with exit code ${GATKRTCCODE}!"
   exit 8
fi

#If RealignerTargetCreator generated output, run IndelRealigner on that output:
if [[ -e ${OUTPUTDIR}${SAMPLE}_realignTargets.list ]]; then
   echo "Running GATK IndelRealigner on ${SAMPLE}_sorted${USEMARKDUP}.bam"
   java -Xmx30g -jar $GATK -T IndelRealigner -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted${USEMARKDUP}.bam -targetIntervals ${OUTPUTDIR}${SAMPLE}_realignTargets.list -o ${OUTPUTDIR}${SAMPLE}${USEMARKDUP}_realigned.bam${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_IR${SAMPLE}.log
   GATKIRCODE=$?
   if [[ $GATKIRCODE -ne 0 ]]; then
      echo "GATK IndelRealigner on ${SAMPLE}_sorted${USEMARKDUP}.bam failed with exit code ${GATKIRCODE}!"
      exit 9
   fi
else
   echo "RealignerTargetCreator failed"
   exit 7
fi
