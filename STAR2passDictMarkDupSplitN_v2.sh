#!/bin/bash
#Read in the command line arguments:
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
if [[ ! -z "$4" ]]; then
   if [[ $4 =~ ^[0-9]+$ ]]; then #If fourth argument exists and is numeric
      NUMPROCS=$4
   else
      READS="${READS} $4" #Append the second read file to the string
      if [[ ! -z "$5" ]]; then
         NUMPROCS=$5
      fi
   fi
fi

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
fi

mkdir -p ${OUTPUTDIR}logs

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Generate index for the reference FASTA:
if [[ ! -e ${REF}.fai ]]; then
   echo "Generating FASTA index for ${REF}"
   $SAMTOOLS faidx ${REF} 2>&1 > ${OUTPUTDIR}logs/samtoolsFaidxRef.log
   FAIDXCODE=$?
   if [[ $FAIDXCODE -ne 0 ]]; then
      echo "samtools faidx on ${REF} failed with exit code ${FAIDXCODE}!"
      exit 2
   fi
else
   echo "Skipping .fai generation for ${REF}"
fi

#Generate sequence dictionary for the reference FASTA (used by GATK):
REFDICT="${REF%.*}.dict"
if [[ ! -e ${REFDICT} ]]; then
   echo "Generating sequence dictionary for ${REF}"
   java -Xmx30g -jar $PICARD CreateSequenceDictionary REFERENCE=${REF} OUTPUT=${REFDICT} 2>&1 > ${OUTPUTDIR}logs/picardDictRef.log
   DICTCODE=$?
   if [[ $DICTCODE -ne 0 ]]; then
      echo "Picard CreateSequenceDictionary on ${REF} failed with exit code ${DICTCODE}!"
      exit 3
   fi
else
   echo "Skipping .dict generation for ${REF}"
fi
#Clear out the STAR temp directory if it exists:
STARTMP="${OUTPUTDIR}${PREFIX}_STARtmp"
if [[ -d "${STARTMP}" ]]; then
   rm -rf ${STARTMP}
fi

#Generate the index for STAR:
GENOMEDIR="${OUTPUTDIR}${PREFIX}_genomeDir"
if [[ -d "${GENOMEDIR}" ]]; then
   rm -rf ${GENOMEDIR}
fi
mkdir -p ${GENOMEDIR}
echo "Generating STAR genome index for ${REF}"
$STAR --runThreadN ${NUMPROCS} --runMode genomeGenerate --genomeDir ${GENOMEDIR} --genomeFastaFiles ${REF} --outFileNamePrefix ${OUTPUTDIR}${PREFIX} --outTmpDir ${STARTMP}
STARIDXCODE=$?
if [[ $STARIDXCODE -ne 0 ]]; then
   echo "STAR genomeGenerate on ${REF} failed with exit code ${STARIDXCODE}!"
   exit 4
fi

#Now run STAR 2-pass (on-the-fly), either with gzipped or uncompressed FASTQ files:
echo "Running STAR 2-pass on ${READS} against ${REF} with ${NUMPROCS} cores"
if [[ ${READS} =~ .gz ]]; then
   $STAR --runThreadN ${NUMPROCS} --genomeDir ${GENOMEDIR} --outTmpDir ${STARTMP} --outSAMtype BAM Unsorted --outStd BAM_Unsorted --outFileNamePrefix ${OUTPUTDIR}${PREFIX}_ --readFilesCommand zcat --readFilesIn ${READS} --twopassMode Basic --outSAMmapqUnique 60 --outSAMattrRGline ID:${PREFIX} LB:${PREFIX} PL:ILLUMINA SM:${PREFIX} 2> ${OUTPUTDIR}logs/STAR2pass_${PREFIX}.stderr | $SAMTOOLS sort -@ ${NUMPROCS} -T ${OUTPUTDIR}${PREFIX}_sorttemp -o ${OUTPUTDIR}${PREFIX}_sorted.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsSort_${PREFIX}.log
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

#Mark duplicates using Picard:
echo "Marking duplicates for ${PREFIX}_sorted.bam"
java -Xmx30g -jar $PICARD MarkDuplicates INPUT=${OUTPUTDIR}${PREFIX}_sorted.bam OUTPUT=${OUTPUTDIR}${PREFIX}_sorted_markdup.bam METRICS_FILE=${OUTPUTDIR}${PREFIX}_markdup_metrics.txt 2>&1 > ${OUTPUTDIR}logs/picardMarkDuplicates${PREFIX}.log
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

#Clean up the mapping results for use by HaplotypeCaller:
echo "Cleaning up mapping results for use by GATK HaplotypeCaller"
echo "using GATK SplitNCigarReads on ${PREFIX}_sorted_markdup.bam"
java -jar $GATK -T SplitNCigarReads -R ${REF} -I ${OUTPUTDIR}${PREFIX}_sorted_markdup.bam -o ${OUTPUTDIR}${PREFIX}_sorted_markdup_splitN.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2>&1 > ${OUTPUTDIR}logs/${PREFIX}_GATK_SplitNCigarReads.log
SPLITNCODE=$?
if [[ $SPLITNCODE -ne 0 ]]; then
   echo "GATK SplitNCigarReads on ${PREFIX}_sorted_markdup.bam failed with exit code ${SPLITNCODE}!"
   exit 8
fi

echo "STAR2passDictMarkDupSplitN on sample ${PREFIX} is complete!"
