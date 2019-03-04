#!/bin/bash
#Read in the command line arguments:
#PREFIX: Prefix to use for all intermediate and output files
# in the pipeline
PREFIX=$1
#REF: Path to the reference used for mapping
REF=$2
#READS: Path to reads for mapping
# For single-end reads, only one file is provided.
# For paired-end reads, both files are provided separated by
# a space, and the whole string is in quotes.
READS=$3
#NUMPROCS: Number of cores to use for mapping
NUMPROCS=8
#SPECIAL: Special options, e.g. interleaved for bwa mem -p
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

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
fi

mkdir -p ${OUTPUTDIR}logs

#Set the paths to the important executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Create a BWA index of the reference FASTA:
if [[ ! -e ${REF}.bwt ]]; then
   echo "Making BWA index of ${REF}"
   $BWA index ${REF} 2>&1 > ${OUTPUTDIR}logs/bwaIndexRef.log
   BWAINDEXCODE=$?
   if [[ $BWAINDEXCODE -ne 0 ]]; then
      echo "BWA indexing of ${REF} failed with exit code ${BWAINDEXCODE}!"
      exit 2
   fi
else
   echo "Skipping BWA index creation for ${REF}"
fi

#Create an index of the reference FASTA:
if [[ ! -e ${REF}.fai ]]; then
   echo "Making FASTA index of ${REF}"
   $SAMTOOLS faidx ${REF} 2>&1 > ${OUTPUTDIR}logs/samtoolsFaidxRef.log
   FAIDXCODE=$?
   if [[ $FAIDXCODE -ne 0 ]]; then
      echo "samtools faidx of ${REF} failed with exit code ${FAIDXCODE}!"
      exit 3
   fi
else
   echo "Skipping .fai creation for ${REF}"
fi

#Create a sequence dictionary of the reference FASTA (for later use by GATK):
REFDICT="${REF%.*}.dict"
if [[ ! -e ${REFDICT} ]]; then
   echo "Creating sequence dictionary for ${REF}"
   java -Xmx30g -jar $PICARD CreateSequenceDictionary REFERENCE=${REF} OUTPUT=${REFDICT} 2>&1 > ${OUTPUTDIR}logs/picardDictRef.log
   DICTCODE=$?
   if [[ $DICTCODE -ne 0 ]]; then
      echo "Picard CreateSequenceDictionary on ${REF} failed with exit code ${DICTCODE}!"
      exit 4
   fi
else
   echo "Skipping .dict creation for ${REF}"
fi

#SPECIAL options may skip alignment and go straight to markdup:
SKIPALN=0
if [[ $SPECIAL =~ "only_markdup" ]]; then #skip to markdup
   SKIPALN=1
fi

#SPECIAL options may add to the extra BWA options:
BWAOPTIONS=""
if [[ $SPECIAL =~ "interleaved" ]]; then #add the -p option
   BWAOPTIONS="${BWAOPTIONS} -p"
fi
if [[ $SPECIAL =~ "comments" ]]; then #add the -C option
   BWAOPTIONS="${BWAOPTIONS} -C"
fi

if [[ "${SKIPALN}" -eq "0" ]]; then
   #Map the reads using BWA mem, convert the output to BAM, and coordinate-sort it using samtools:
   echo "Mapping ${READS} to ${REF} with bwa mem using ${NUMPROCS} cores with extra options ${BWAOPTIONS}"
   $BWA mem ${BWAOPTIONS} -t ${NUMPROCS} -R "@RG\\tID:${PREFIX}\\tSM:${PREFIX}\\tPL:ILLUMINA\\tLB:${PREFIX}" ${REF} ${READS} 2> ${OUTPUTDIR}logs/bwaMem${PREFIX}.stderr | $SAMTOOLS view -u -@ ${NUMPROCS} - 2> ${OUTPUTDIR}logs/samtoolsView${PREFIX}.stderr | $SAMTOOLS sort -@ ${NUMPROCS} -m 3G -T ${OUTPUTDIR}${PREFIX}_sorttemp -o ${OUTPUTDIR}${PREFIX}_sorted.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsSort${PREFIX}.log
   BWAMEMCODE=$?
   if [[ $BWAMEMCODE -ne 0 ]]; then
      echo "bwa mem of ${READS} against ${REF} failed with exit code ${BWAMEMCODE}!"
      exit 5
   fi
else
   echo "Skipping alignment as requested by ${SPECIAL} for sample ${PREFIX}"
fi

#Mark duplicates using Picard:
echo "Marking duplicates using Picard for sample ${PREFIX}"
java -Xmx30g -jar $PICARD MarkDuplicates INPUT=${OUTPUTDIR}${PREFIX}_sorted.bam OUTPUT=${OUTPUTDIR}${PREFIX}_sorted_markdup.bam METRICS_FILE=${OUTPUTDIR}${PREFIX}_markdup_metrics.txt 2>&1 > ${OUTPUTDIR}logs/picardMarkDuplicates${PREFIX}.log
MARKDUPCODE=$?
if [[ $MARKDUPCODE -ne 0 ]]; then
   echo "Picard MarkDuplicates on ${PREFIX}_sorted.bam failed with exit code ${MARKDUPCODE}!"
   exit 6
fi

#Index the BAM produced by Picard MarkDuplicates:
echo "Indexing the duplicate-marked BAM for sample ${PREFIX}"
$SAMTOOLS index ${OUTPUTDIR}${PREFIX}_sorted_markdup.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsIndex${PREFIX}.log
BAMINDEXCODE=$?
if [[ $BAMINDEXCODE -ne 0 ]]; then
   echo "samtools index on ${PREFIX}_sorted_markdup.bam failed with exit code ${BAMINDEXCODE}!"
   exit 7
fi

echo "iADMD complete for sample ${PREFIX}!"
