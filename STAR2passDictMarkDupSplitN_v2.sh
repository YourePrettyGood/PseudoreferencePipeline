#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [number of cores] [special options]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and .dict.\n"
   printf "Number of cores allocated is 8 by default (match with --cpus-per-task).\n"
   printf "Special options is a comma-separated list of flags to deviate\n"
   printf " from default operation.\n"
   printf "'no_markdup' skips the Picard MarkDuplicates step.\n"
   printf "'no_splitN' skips the GATK SplitNCigarReads step.\n"
   printf "'intronMotif' adds the XS tag for reads aligning across canonical\n"
   printf " splice junctions (nominally for compatibility with Cufflinks and\n"
   printf " StringTie).\n"
   printf "'mem_#[mg]' specifies the maximum memory allowed for Java during\n"
   printf " Picard MarkDuplicates and SplitNCigarReads.\n"
   printf " (default is 30g for 30 GB).\n"
   printf " The #[mg] is appended to -Xmx in the java call.\n"
   printf "'cleanup' omits primary analysis functions and instead deletes all non-log intermediate files.\n"
   exit 1
fi

#PREFIX: Prefix used for all intermediate and output files in the
# pipeline
PREFIX=$1
#REF: The reference FASTA used for mapping
REF=$2
#READS: Path to the reads files used for mapping
# For single-end reads, only one file is provided.
# For paired-end reads, the two files are separated by a space,
# and the whole string is in quotes.
READS=$3
#NUMPROCS: Number of cores to use for mapping
NUMPROCS=8
#SPECIAL: Special options, e.g.
SPECIAL=""
if [[ ! -z "$4" ]]; then
   if [[ $4 =~ ^[0-9]+$ ]]; then #If fourth argument exists and is numeric
      NUMPROCS=$4
      if [[ ! -z "$5" ]]; then #If fifth argument exists, it is SPECIAL here
         SPECIAL="$5"
      fi
   else
      READS="${READS} $4" #Append the second read file to the string
      if [[ ! -z "$5" ]]; then
         NUMPROCS=$5
         if [[ ! -z "$6" ]]; then #If sixth argument exists, it is SPECIAL here
            SPECIAL="$6"
         fi
      fi
   fi
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
echo "Running java steps with -Xmx${JAVAMEM} for sample ${PREFIX}"

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
fi

mkdir -p ${OUTPUTDIR}logs

#SPECIAL options may indicate cleanup of intermediate files,
# but not log files:
if [[ $SPECIAL =~ "cleanup" ]]; then
   rm -f ${OUTPUTDIR}${PREFIX}_sorted.bam
   if [[ ! $SPECIAL =~ "no_markdup" ]]; then
      OUTPUTBAM="${OUTPUTDIR}${PREFIX}_sorted_markdup_splitN.bam"
      rm -f ${OUTPUTDIR}${PREFIX}_sorted_markdup.bam ${OUTPUTDIR}${PREFIX}_sorted_markdup.bam.bai
   else
      OUTPUTBAM="${OUTPUTDIR}${PREFIX}_sorted_splitN.bam"
   fi
   if [[ ! $SPECIAL =~ "no_splitN" ]]; then
      rm -f ${OUTPUTBAM} ${OUTPUTBAM}.bai
   fi
   echo "Cleanup complete for sample ${PREFIX}"
   exit 0
fi

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Check for FASTA index of the reference FASTA:
if [[ ! -e ${REF}.fai ]]; then
   echo "Missing FASTA index of ${REF}, please run indexDictFai.sh"
   exit 2
fi

#Check for sequence dictionary of the reference FASTA (used by GATK):
REFDICT="${REF%.*}.dict"
if [[ ! -e ${REFDICT} ]]; then
   echo "Missing sequence dictionary for ${REF}, please run indexDictFai.sh"
   exit 3
fi

#Clear out the STAR temp directory if it exists:
STARTMP="${OUTPUTDIR}${PREFIX}_STARtmp"
if [[ -d "${STARTMP}" ]]; then
   rm -rf ${STARTMP}
fi

#Check for the index for STAR:
GENOMEDIR="${REF}_genomeDir"
if [[ ! -d "${GENOMEDIR}" ]]; then
   echo "Missing STAR index of ${REF}, please run indexDictFai.sh STAR"
   exit 4
fi

#Deal with SPECIAL for STAR options, like intronMotif:
STAROPTIONS=""
if [[ "${SPECIAL}" =~ "intronMotif" ]]; then
   STAROPTIONS="${STAROPTIONS} --outSAMstrandField intronMotif"
fi

#Check for necessary scripts/programs:
if [[ ! -x "$(command -v ${SAMTOOLS})" ]]; then
   echo "SAMtools appears to be missing, could not find at SAMTOOLS=${SAMTOOLS}."
   exit 18;
fi
if [[ ! -e ${PICARD} ]]; then
   echo "Picard appears to be missing, could not find at PICARD=${PICARD}."
   exit 11;
fi

#Now run STAR 2-pass (on-the-fly), either with gzipped or uncompressed FASTQ files:
echo "Running STAR 2-pass on ${READS} against ${REF} with ${NUMPROCS} cores"
if [[ ${READS} =~ ".gz" ]]; then
   $STAR --runThreadN ${NUMPROCS} --genomeDir ${GENOMEDIR} --outTmpDir ${STARTMP} --outSAMtype BAM Unsorted --outStd BAM_Unsorted --outFileNamePrefix ${OUTPUTDIR}${PREFIX}_ --readFilesCommand zcat --readFilesIn ${READS} --twopassMode Basic --outSAMmapqUnique 60 ${STAROPTIONS} --outSAMattrRGline ID:${PREFIX} LB:${PREFIX} PL:ILLUMINA SM:${PREFIX} 2> ${OUTPUTDIR}logs/STAR2pass_${PREFIX}.stderr | $SAMTOOLS sort -@ ${NUMPROCS} -T ${OUTPUTDIR}${PREFIX}_sorttemp -o ${OUTPUTDIR}${PREFIX}_sorted.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsSort_${PREFIX}.log
   STARCODE=$?
   if [[ $STARCODE -ne 0 ]]; then
      echo "STAR mapping of ${READS} on ${REF} failed with exit code ${STARCODE}!"
      exit 5
   fi
else
   $STAR --runThreadN ${NUMPROCS} --genomeDir ${GENOMEDIR} --outTmpDir ${STARTMP} --outSAMtype BAM Unsorted --outStd BAM_Unsorted --outFileNamePrefix ${OUTPUTDIR}${PREFIX}_ --readFilesIn ${READS} --twopassMode Basic --outSAMmapqUnique 60 --outSAMattrRGline ID:${PREFIX} LB:${PREFIX} PL:ILLUMINA SM:${PREFIX} 2> ${OUTPUTDIR}logs/STAR2pass_${PREFIX}.stderr | $SAMTOOLS sort -@ ${NUMPROCS} -T ${OUTPUTDIR}${PREFIX}_sorttemp -o ${OUTPUTDIR}${PREFIX}_sorted.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsSort_${PREFIX}.log
   STARCODE=$?
   if [[ $STARCODE -ne 0 ]]; then
      echo "STAR mapping of ${READS} on ${REF} failed with exit code ${STARCODE}!"
      exit 5
   fi
fi

if [[ ${SPECIAL} =~ "no_markdup" ]]; then
   echo "Skipping marking of duplicates for sample ${PREFIX} due to special flag ${SPECIAL}"
else
   #Mark duplicates using Picard:
   echo "Marking duplicates for ${PREFIX}_sorted.bam"
   java -Xmx${JAVAMEM} -jar $PICARD MarkDuplicates INPUT=${OUTPUTDIR}${PREFIX}_sorted.bam OUTPUT=${OUTPUTDIR}${PREFIX}_sorted_markdup.bam METRICS_FILE=${OUTPUTDIR}${PREFIX}_markdup_metrics.txt 2>&1 > ${OUTPUTDIR}logs/picardMarkDuplicates${PREFIX}.log
   MARKDUPCODE=$?
   if [[ $MARKDUPCODE -ne 0 ]]; then
      echo "Picard MarkDuplicates on ${PREFIX}_sorted.bam failed with exit code ${MARKDUPCODE}!"
      exit 6
   fi
   
   #Index the BAM produced by Picard:
   echo "Indexing the duplicate-marked BAM ${PREFIX}_sorted_markdup.bam"
   $SAMTOOLS index ${OUTPUTDIR}${PREFIX}_sorted_markdup.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsIndex${PREFIX}.log
   BAMIDXCODE=$?
   if [[ $BAMIDXCODE -ne 0 ]]; then
      echo "samtools index on ${PREFIX}_sorted_markdup.bam failed with exit code ${BAMIDXCODE}!"
      exit 7
   fi
fi

if [[ ${SPECIAL} =~ "no_splitN" ]]; then
   echo "Skipping GATK SplitNCigarReads for sample ${PREFIX} due to special options ${SPECIAL}"
else
   if [[ ! -e "${GATK}" ]]; then
      echo "GATK appears to be missing, could not find at GATK=${GATK}."
      exit 20;
   fi
   INPUTBAM="${OUTPUTDIR}${PREFIX}_sorted_markdup.bam"
   OUTPUTBAM="${OUTPUTDIR}${PREFIX}_sorted_markdup_splitN.bam"
   if [[ ${SPECIAL} =~ "no_markdup" ]]; then
      INPUTBAM="${OUTPUTDIR}${PREFIX}_sorted.bam"
      OUTPUTBAM="${OUTPUTDIR}${PREFIX}_sorted_splitN.bam"
      ${SAMTOOLS} index ${INPUTBAM}
      NOMDINDEXCODE=$?
      if [[ $NOMDINDEXCODE -ne 0 ]]; then
         echo "Indexing of non-markdup BAM before SplitNCigarReads failed with exit code ${NOMDINDEXCODE} for sample ${PREFIX}"
         exit 9
      fi
   fi
   #Clean up the mapping results for use by HaplotypeCaller:
   echo "Cleaning up mapping results for use by GATK HaplotypeCaller"
   echo "using GATK SplitNCigarReads on ${INPUTBAM}"
   java -Xmx${JAVAMEM} -jar $GATK -T SplitNCigarReads -R ${REF} -I ${INPUTBAM} -o ${OUTPUTBAM} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2>&1 > ${OUTPUTDIR}logs/${PREFIX}_GATK_SplitNCigarReads.log
   SPLITNCODE=$?
   if [[ $SPLITNCODE -ne 0 ]]; then
      echo "GATK SplitNCigarReads on ${PREFIX}_sorted_markdup.bam failed with exit code ${SPLITNCODE}!"
      exit 8
   fi
fi

echo "STAR2passDictMarkDupSplitN on sample ${PREFIX} is complete!"
