#!/bin/bash
#Read in the command line arguments:
#Output usage if none supplied:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [extra option]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Extra option is either 'misencoded' or 'filter_bases_not_stored'\n"
   printf "where 'misencoded' forces GATK to convert PHRED+64 qual scores to\n"
   printf "PHRED+33 qual scores, and 'filter_bases_not_stored' forces GATK to\n"
   printf "ignore alignments where the seq or qual columns are *.\n"
   exit 1
fi
#SAMPLE: Prefix used for all intermediate and output files of the pipeline
SAMPLE=$1
#REFERENCE: Path to the FASTA of the reference used for mapping
REFERENCE=$2
#MISENCODED: Flag used to trigger conversion of PHRED+64 quality scores in
# your BAM to PHRED+33, because GATK does not auto-convert PHRED+64, and
# complains if you give it PHRED+64 data...
#FILTERZEROLEN: Flag used to trigger ignoring of alignments where bases are
# not stored in the BAM (potentially for zero-length reads?)
if [[ -z "$3" ]]; then
   MISENCODED=""
   FILTERZEROLEN=""
elif [[ $3 =~ "filter_bases_not_stored" ]]; then
   FILTERZEROLEN=" --filter_bases_not_stored"
   MISENCODED=""
elif [[ $3 =~ "filter_mismatching_base_and_quals" ]]; then
   FILTERZEROLEN=" --filter_mismatching_base_and_quals"
   MISENCODED=""
else
   FILTERZEROLEN=""
   MISENCODED=" --fix_misencoded_quality_scores"
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

if [[ ! -e ${GATK} ]]; then
   echo "Path to GATK is invalid."
   exit 2
fi
if [[ ! -e ${SAMTOOLS} ]]; then
   echo "Path to Samtools is invalid."
   exit 3
fi

#Verify that input files exist:
if [[ ! -e ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam ]]; then
   echo "Input BAM file ${SAMPLE}_sorted_markdup_splitN.bam is missing.  Cannot proceed."
   exit 4
fi
if [[ ! -e ${REFERENCE} ]]; then
   echo "Reference FASTA ${REFERENCE} is missing.  Cannot proceed."
   exit 5
fi

#Check that the input BAM has an index, since GATK complains if it's missing:
if [[ ! -e ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam.bai ]]; then
   echo "Index was missing from input BAM.  Now indexing."
   echo "$SAMTOOLS index ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam"
   $SAMTOOLS index ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam
   BAMIDXCODE=$?
   if [[ $BAMIDXCODE -ne 0 ]]; then
      echo "samtools index on ${SAMPLE}_sorted_markdup_splitN.bam failed with exit code ${BAMIDXCODE}!"
      exit 6
   fi
fi

#Run GATK's RealignerTargetCreator to generate a list of targets for realignment:
echo "Running RealignerTargetCreator for sample ${SAMPLE}:"
echo "java -Xmx30g -jar $GATK -T RealignerTargetCreator -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam -o ${OUTPUTDIR}${SAMPLE}_realignTargets.list${MISENCODED}${FILTERZEROLEN} 2>&1 > ${OUTPUTDIR}logs/GATK_RTC${SAMPLE}.log"
java -Xmx30g -jar $GATK -T RealignerTargetCreator -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam -o ${OUTPUTDIR}${SAMPLE}_realignTargets.list${MISENCODED}${FILTERZEROLEN} 2>&1 > ${OUTPUTDIR}logs/GATK_RTC${SAMPLE}.log
GATKRTCCODE=$?
if [[ $GATKRTCCODE -ne 0 ]]; then
   echo "GATK RealignerTargetCreator on ${SAMPLE}_sorted_markdup_splitN.bam failed with exit code ${GATKRTCCODE}!"
   exit 7
fi

#If RealignerTargetCreator generated output, run IndelRealigner on that output:
echo "Running IndelRealigner on ${SAMPLE}_sorted_markdup_splitN.bam:"
if [[ -e ${OUTPUTDIR}${SAMPLE}_realignTargets.list ]]; then
   echo "java -Xmx30g -jar $GATK -T IndelRealigner -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam -targetIntervals ${OUTPUTDIR}${SAMPLE}_realignTargets.list -o ${OUTPUTDIR}${SAMPLE}_markdup_realigned.bam${FILTERZEROLEN}${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_IR${SAMPLE}.log"
   java -Xmx30g -jar $GATK -T IndelRealigner -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam -targetIntervals ${OUTPUTDIR}${SAMPLE}_realignTargets.list -o ${OUTPUTDIR}${SAMPLE}_markdup_realigned.bam${FILTERZEROLEN}${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_IR${SAMPLE}.log
   GATKIRCODE=$?
   if [[ $GATKIRCODE -ne 0 ]]; then
      echo "GATK IndelRealigner on ${SAMPLE}_sorted_markdup_splitN.bam failed with exit code ${GATKIRCODE}!"
      exit 8
   fi
else
   echo "RealignerTargetCreator failed"
   exit 7
fi
