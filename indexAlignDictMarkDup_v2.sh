#!/bin/bash
#Read in the command line arguments:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <prefix> <reference genome> [read file(s)] [number of cores] [special options]\n"
   printf "Prefix is used to standardize intermediate and output file names.\n"
   printf "Reference genome must be wrapped FASTA with .fai and BWA index.\n"
   printf "Number of cores allocated is 8 by default (match with --cpus-per-task).\n"
   printf "Special options is a comma-separated list of flags to deviate\n"
   printf " from default operation.\n"
   printf "'only_bwa' omits the Picard MarkDuplicates step.\n"
   printf " (Make sure to use 'no_markdup' for any later steps if you do this.)\n"
   printf "'no_markdup' does the same exact thing as 'only_bwa'\n"
   printf "'only_markdup' omits the mapping step and only performs marking of\n"
   printf " marking of duplicates with Picard MarkDuplicates.\n"
   printf "'mem_#[mg]' specifies the maximum memory allowed for Java during\n"
   printf " Picard MarkDuplicates (default is 30g for 30 GB).\n"
   printf " The #[mg] is appended to -Xmx in the java call.\n"
   printf "'interleaved' tells BWA to interpret the single read file as an\n"
   printf " interleaved FASTQ, i.e. odd-indexed reads are R1, even-indexed\n"
   printf " reads are R2 of a paired-end run.\n"
   printf " This can be useful for processing e.g. parsed reads from 10x\n"
   printf " LongRanger Basic (along with 'comments' to preserve 10x barcodes).\n"
   printf "'comments' tells BWA to include string parts of the FASTQ header\n"
   printf " after the first space in the output BAM as a comment to the alignment.\n"
   printf " This is useful to preserve 10x barcodes from a FASTQ produced by\n"
   printf " 10x LongRanger Basic.\n"
   printf "'all_alignments' indicates using the '-a' option for BWA, outputting\n"
   printf " all found alignments for a given read (extra alignments are marked\n"
   printf " as secondary).\n"
   printf "'cleanup' omits primary analysis functions and instead deletes all non-log intermediate files.\n"
   exit 1
fi

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

#Notes for future addition of Bowtie2:
#Add option to specify BT2 instead of BWA-MEM
#Default options (maybe?):
#--local --sensitive-local (or --very-sensitive-local) -t (for logging)
#Consider adding --xeq for extended CIGAR support
#SPECIAL =~ 'interleaved' => BT2READS="--interleaved ${READS}"
#SPECIAL =~ 'all_alignments' => BT2OPTIONS="${BT2OPTIONS} -a"
#If ${READS} has one token => BT2READS="-U ${READS}"
#Else split into tokens => BT2READS="-1 ${readsarr[0]} -2 ${readsarr[1]}"
#bowtie2 ${BT2OPTIONS} -p ${NUMPROCS} -rg-id "${PREFIX}" -rg "SM:${PREFIX}" -rg "PL:ILLUMINA" -rg "LB:${PREFIX}" -x ${REF} ${BT2READS} 2> ${OUTPUTDIR}logs/bt2${PREFIX}.stderr | $SAMTOOLS view -u -@ ${NUMPROCS} - 2> ${OUTPUTDIR}logs/samtoolsView${PREFIX}.stderr | $SAMTOOLS sort -@ ${NUMPROCS} -m 3G -T ${OUTPUTDIR}${PREFIX}_sorttemp -o ${OUTPUTDIR}${PREFIX}_sorted.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsSort${PREFIX}.log

#Notes for future addition of minimap2:
#Add option to specify minimap2
#Default options:
#-ax sr
#Maybe --cs or --cs=long? Not necessary for typical variant callers though
#Maybe --eqx to support extended CIGAR
#SPECIAL =~ 'interleaved' => No options necessary
#SPECIAL =~ 'all_alignments' => maybe -p 0.0? not sure
#SPECIAL =~ 'comments' => MM2OPTIONS="${MM2OPTIONS} -y"
#minimap2 ${MM2OPTIONS} -t ${NUMPROCS} -R "@RG\\tID:${PREFIX}\\tSM:${PREFIX}\\tPL:ILLUMINA\\tLB:${PREFIX}" ${REF} ${READS} 2> ${OUTPUTDIR}logs/minimap2${PREFIX}.stderr | $SAMTOOLS view -u -@ ${NUMPROCS} - 2> ${OUTPUTDIR}logs/samtoolsView${PREFIX}.stderr | $SAMTOOLS sort -@ ${NUMPROCS} -m 3G -T ${OUTPUTDIR}${PREFIX}_sorttemp -o ${OUTPUTDIR}${PREFIX}_sorted.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsSort${PREFIX}.log

#Check if memory for GATK/java is specified:
#Format of SPECIAL argument is mem_[0-9]+[MGmg]? 
#The part after the underscore gets used as a suffix to -Xmx
#Since we're using a regex to capture, this is not an attack vector
JAVAMEM="30g"
JAVAMEMRE='(mem)_([0-9]+[MGmg]?)'
if [[ ${SPECIAL} =~ $JAVAMEMRE ]]; then
   JAVAMEM="${BASH_REMATCH[2]}"
fi
echo "Running java steps with -Xmx${JAVAMEM} for sample ${SAMPLE}"

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
fi

mkdir -p ${OUTPUTDIR}logs

#Set the paths to the important executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#SPECIAL options may skip alignment and go straight to markdup,
# or perform only alignment and skip markdup:
SKIPALN=0
SKIPMD=0
if [[ $SPECIAL =~ "only_markdup" ]]; then #skip to markdup
   SKIPALN=1
fi
if [[ $SPECIAL =~ "only_bwa" || $SPECIAL =~ "no_markdup" ]]; then #skip markdup
   SKIPMD=1
fi

#SPECIAL options may indicate cleanup of intermediate files,
# but not log files:
if [[ $SPECIAL =~ "cleanup" ]]; then
   if [[ "${SKIPALN}" -eq "0" ]]; then
      rm -f ${OUTPUTDIR}${PREFIX}_sorted.bam ${OUTPUTDIR}${PREFIX}_sorted.bam.bai
   fi
   if [[ "${SKIPMD}" -eq "0" ]]; then
      rm -f ${OUTPUTDIR}${PREFIX}_sorted_markdup.bam ${OUTPUTDIR}${PREFIX}_sorted_markdup.bam.bai
   fi
   echo "Cleanup complete for sample ${PREFIX}"
   exit 0
fi

#Check for a BWA index of the reference FASTA:
if [[ ! -e ${REF}.bwt ]]; then
   echo "Missing BWA index of ${REF}, please run indexDictFai.sh"
   exit 2
fi

#Check for an index of the reference FASTA:
if [[ ! -e ${REF}.fai ]]; then
   echo "Missing FASTA index of ${REF}, please run indexDictFai.sh"
   exit 3
fi

#Check for a sequence dictionary of the reference FASTA (for later use by GATK):
REFDICT="${REF%.*}.dict"
if [[ ! -e ${REFDICT} ]]; then
   echo "Missing sequence dictionary for ${REF}, please run indexDictFai.sh"
   exit 4
fi

#SPECIAL options may add to the extra BWA options:
BWAOPTIONS=""
if [[ $SPECIAL =~ "interleaved" ]]; then #add the -p option
   BWAOPTIONS="${BWAOPTIONS} -p"
fi
if [[ $SPECIAL =~ "comments" ]]; then #add the -C option
   BWAOPTIONS="${BWAOPTIONS} -C"
fi
if [[ $SPECIAL =~ "all_alignments" ]]; then #add the -a option
   BWAOPTIONS="${BWAOPTIONS} -a"
fi

if [[ ! -x "$(command -v ${BWA})" ]]; then
   echo "BWA appears to be missing, could not find at BWA=${BWA}."
   exit 16;
fi
if [[ ! -x "$(command -v ${SAMTOOLS})" ]]; then
   echo "SAMtools appears to be missing, could not find at SAMTOOLS=${SAMTOOLS}."
   exit 18;
fi
if [[ ! -e ${PICARD} ]]; then
   echo "Picard appears to be missing, could not find at PICARD=${PICARD}."
   exit 11;
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

#Index the BAM produced by BWA+SAMtools:
echo "Indexing the sorted BAM for sample ${PREFIX}"
$SAMTOOLS index ${OUTPUTDIR}${PREFIX}_sorted.bam 2>&1 > ${OUTPUTDIR}logs/samtoolsIndexSorted${PREFIX}.log
INDEXCODE=$?
if [[ $INDEXCODE -ne 0 ]]; then
   echo "samtools index on ${PREFIX}_sorted.bam failed with exit code ${INDEXCODE}!"
   exit 8
fi


if [[ "${SKIPMD}" -eq "0" ]]; then
   #Mark duplicates using Picard:
   echo "Marking duplicates using Picard for sample ${PREFIX}"
   java -Xmx${JAVAMEM} -jar $PICARD MarkDuplicates INPUT=${OUTPUTDIR}${PREFIX}_sorted.bam OUTPUT=${OUTPUTDIR}${PREFIX}_sorted_markdup.bam METRICS_FILE=${OUTPUTDIR}${PREFIX}_markdup_metrics.txt 2>&1 > ${OUTPUTDIR}logs/picardMarkDuplicates${PREFIX}.log
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
else
   echo "Skipping markdup as requested by ${SPECIAL} for sample ${PREFIX}"
fi

echo "iADMD complete for sample ${PREFIX}!"
