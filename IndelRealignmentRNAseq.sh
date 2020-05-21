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
   printf "'filter_bases_not_stored' forces GATK to skip alignments where\n"
   printf " the seq or qual columns are *.\n"
   printf "'filter_mismatching_base_and_quals' forces GATK to skip alignments\n"
   printf " where the seq and qual columns have different lengths.\n"
   printf "These latter two are generally unnecessary, as they indicate a\n"
   printf " malformatted or truncated BAM (thus upstream diagnosis is needed).\n"
   printf "'mem_#[mg]' specifies the maximum memory allowed for Java\n"
   printf " (default is 30g for 30 GB).\n"
   printf " The #[mg] is appended to -Xmx in the java call.\n"
   printf "'cleanup' omits primary analysis functions and instead deletes all non-log intermediate files.\n"
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
SPECIAL=$3
MISENCODED=""
FILTERZEROLEN=""
if [[ ${SPECIAL} =~ "filter_bases_not_stored" ]]; then
   FILTERZEROLEN=" --filter_bases_not_stored"
fi
if [[ ${SPECIAL} =~ "filter_mismatching_base_and_quals" ]]; then
   FILTERZEROLEN=" --filter_mismatching_base_and_quals"
fi
if [[ ${SPECIAL} =~ "misencoded" ]]; then
   MISENCODED=" --fix_misencoded_quality_scores"
fi
#no_markdup: Specify using a BAM that hasn't been markduped
#if [[ ${SPECIAL} =~ "no_markdup" ]]; then
#   MARKDUP=""
#   NOMARKDUP="_nomarkdup"
#else
#   MARKDUP="_markdup"
#   NOMARKDUP=""
#fi

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
   rm -f ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam.bai ${OUTPUTDIR}${SAMPLE}_realignTargets.list ${OUTPUTDIR}${SAMPLE}_markdup_realigned.bam ${OUTPUTDIR}${SAMPLE}_markdup_realigned.bai
   echo "Cleanup complete for sample ${SAMPLE}"
   exit 0
fi

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

if [[ ! -e "${GATK}" ]]; then
   echo "GATK appears to be missing, could not find at GATK=${GATK}."
   exit 20;
fi
if [[ ! -x "$(command -v ${SAMTOOLS})" ]]; then
   echo "SAMtools appears to be missing, could not find at SAMTOOLS=${SAMTOOLS}."
   exit 18;
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
echo "java -Xmx${JAVAMEM} -jar $GATK -T RealignerTargetCreator -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam -o ${OUTPUTDIR}${SAMPLE}_realignTargets.list${MISENCODED}${FILTERZEROLEN} 2>&1 > ${OUTPUTDIR}logs/GATK_RTC${SAMPLE}.log"
java -Xmx${JAVAMEM} -jar $GATK -T RealignerTargetCreator -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam -o ${OUTPUTDIR}${SAMPLE}_realignTargets.list${MISENCODED}${FILTERZEROLEN} 2>&1 > ${OUTPUTDIR}logs/GATK_RTC${SAMPLE}.log
GATKRTCCODE=$?
if [[ $GATKRTCCODE -ne 0 ]]; then
   echo "GATK RealignerTargetCreator on ${SAMPLE}_sorted_markdup_splitN.bam failed with exit code ${GATKRTCCODE}!"
   exit 7
fi

#If RealignerTargetCreator generated output, run IndelRealigner on that output:
echo "Running IndelRealigner on ${SAMPLE}_sorted_markdup_splitN.bam:"
if [[ -e ${OUTPUTDIR}${SAMPLE}_realignTargets.list ]]; then
   echo "java -Xmx${JAVAMEM} -jar $GATK -T IndelRealigner -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam -targetIntervals ${OUTPUTDIR}${SAMPLE}_realignTargets.list -o ${OUTPUTDIR}${SAMPLE}_markdup_realigned.bam${FILTERZEROLEN}${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_IR${SAMPLE}.log"
   java -Xmx${JAVAMEM} -jar $GATK -T IndelRealigner -R ${REFERENCE} -I ${OUTPUTDIR}${SAMPLE}_sorted_markdup_splitN.bam -targetIntervals ${OUTPUTDIR}${SAMPLE}_realignTargets.list -o ${OUTPUTDIR}${SAMPLE}_markdup_realigned.bam${FILTERZEROLEN}${MISENCODED} 2>&1 > ${OUTPUTDIR}logs/GATK_IR${SAMPLE}.log
   GATKIRCODE=$?
   if [[ $GATKIRCODE -ne 0 ]]; then
      echo "GATK IndelRealigner on ${SAMPLE}_sorted_markdup_splitN.bam failed with exit code ${GATKIRCODE}!"
      exit 8
   fi
else
   echo "RealignerTargetCreator failed"
   exit 7
fi
