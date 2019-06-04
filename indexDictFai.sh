#!/bin/bash
#Read in the command line arguments:
#REF: Path to the reference used for mapping
REF=$1

mkdir -p logs

#Set the paths to the important executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Check that the reference FASTA is wrapped (heuristic):
if [[ `head -n2 ${REF} | tail -n1 | tr -d "\n" | wc -c` -gt 1000 ]]; then
   echo "It appears that ${REF} is not wrapped.  Please wrap it, or else GATK will complain."
   exit 1;
fi

#Create a BWA index of the reference FASTA:
if [[ ! -e ${REF}.bwt ]]; then
   echo "Making BWA index of ${REF}"
   $BWA index ${REF} 2>&1 > logs/bwaIndexRef.log
else
   echo "Skipping BWA index creation for ${REF}"
fi

#Make STAR index if requested:
if [[ "$#" -gt 1 ]]; then
   if [[ ! -d "${REF}_genomeDir" ]]; then
      echo "Making STAR index of ${REF}"
      mkdir -p ${REF}_genomeDir
      $STAR --runThreadN 1 --runMode genomeGenerate --genomeDir ${REF}_genomeDir --genomeFastaFiles ${REF}
      STARIDXCODE=$?
      if [[ $STARIDXCODE -ne 0 ]]; then
         echo "STAR genomeGenerate on ${REF} failed with exit code ${STARIDXCODE}!"
         exit 4
      fi
   else
      echo "Skipping STAR index creation for ${REF}"
   fi
fi


#Create an index of the reference FASTA:
if [[ ! -e ${REF}.fai ]]; then
   echo "Making FASTA index of ${REF}"
   $SAMTOOLS faidx ${REF} 2>&1 > logs/samtoolsFaidxRef.log
else
   echo "Skipping .fai creation for ${REF}"
fi

#Create a sequence dictionary of the reference FASTA (for later use by GATK):
REFDICT="${REF%.*}.dict"
if [[ ! -e ${REFDICT} ]]; then
   echo "Creating sequence dictionary for ${REF}"
   java -Xmx30g -jar $PICARD CreateSequenceDictionary REFERENCE=${REF} OUTPUT=${REFDICT} 2>&1 > logs/picardDictRef.log
else
   echo "Skipping .dict creation for ${REF}"
fi

if [[ "$#" -gt 1 ]]; then
   echo "Done creating BWA index, STAR index, .dict, and .fai for ${REF}"
else
   echo "Done creating BWA index, .dict, and .fai for ${REF}"
fi
