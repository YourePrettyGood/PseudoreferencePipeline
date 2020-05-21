#!/bin/bash
#Read in the command line arguments:
#REF: Path to the reference used for mapping
REF=$1
#MAPPER: Which mapper to prepare reference for
# i.e. BWA (default), STAR, BT2, MM2, HISAT2 (more to be added later)
MAPPER=$2
#ANNOTATION: Optional, annotation to supply to STAR
ANNOTATION=$3

#JAVAMEM: Hard-coded for now, memory to allocate for Picard:
JAVAMEM="-Xmx4g"

mkdir -p logs

#Set the paths to the important executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Check that the reference FASTA is wrapped (heuristic):
if [[ `head -n2 ${REF} | tail -n1 | tr -d "\n" | wc -c` -gt 1000 ]]; then
   echo "It appears that ${REF} is not wrapped.  Please wrap it, or else GATK will complain."
   exit 1;
fi

if [[ "$#" -le "1" ]]; then
   MAPPER="BWA";
fi

#Notes for future addition of Bowtie2:
#Add some way to specify BT2 vs. BWA vs. STAR
#bowtie2-build -f ${REF} ${REF} 2>&1 > logs/bt2IndexRef.log

#Notes for future addition of minimap2:
#As above about specifying which mapper
#minimap2 -d ${REF}.mmi ${REF} 2>&1 > logs/minimap2IndexRef.log

if [[ "${MAPPER}" == "BWA" ]]; then
   if [[ ! -x "$(command -v ${BWA})" ]]; then
      echo "BWA appears to be missing, could not find at BWA=${BWA}."
      exit 16;
   fi
   #Create a BWA index of the reference FASTA:
   if [[ ! -e ${REF}.bwt ]]; then
      echo "Making BWA index of ${REF}"
      $BWA index ${REF} 2>&1 > logs/bwaIndexRef.log
      BWAIDXCODE=$?
      if [[ $BWAIDXCODE -ne 0 ]]; then
         echo "bwa index on ${REF} failed with exit code ${BWAIDXCODE}!"
         exit 4
      fi
   else
      echo "Skipping BWA index creation for ${REF}"
   fi
elif [[ "${MAPPER}" == "STAR" ]]; then
   if [[ ! -x "$(command -v ${STAR})" ]]; then
      echo "STAR appears to be missing, could not find at STAR=${STAR}."
      exit 15;
   fi
   #Make STAR index if requested:
   if [[ ! -d "${REF}_genomeDir" ]]; then
      ANNOTATIONARGS=""
      if [[ ! -z "${ANNOTATION}" ]]; then
         ANNOTATIONARGS="--sjdbGTFfile ${ANNOTATION}"
         #Assume GTF by default, but GFF if extension is .gff or .gff3:
         if [[ "${ANNOTATION}" =~ "\.gff3?$" ]]; then
            ANNOTATIONARGS="${ANNOTATIONARGS} --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon exon"
         else
            ANNOTATIONARGS="${ANNOTATIONARGS} --sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFfeatureExon exon"
         fi
      fi
      echo "Making STAR index of ${REF}"
      mkdir -p ${REF}_genomeDir
      ${STAR} --runThreadN 1 --runMode genomeGenerate --genomeDir ${REF}_genomeDir --genomeFastaFiles ${REF} ${ANNOTATIONARGS} 2>&1 > logs/starGenomeGenerate.log
      STARIDXCODE=$?
      if [[ $STARIDXCODE -ne 0 ]]; then
         echo "STAR genomeGenerate on ${REF} failed with exit code ${STARIDXCODE}!"
         exit 4
      fi
   else
      echo "Skipping STAR index creation for ${REF}"
   fi
elif [[ "${MAPPER}" == "BT2" ]]; then
   if [[ ! -x "$(command -v ${BT2BUILD})" ]]; then
      echo "bowtie2-build appears to be missing, could not find at BT2BUILD=${BT2BUILD}."
      exit 14;
   fi
   #Create a bowtie2 index of the reference FASTA:
   if [[ ! -e ${REF}.1.bt2 ]]; then
      echo "Making bowtie2 index of ${REF}"
      ${BT2BUILD} -f ${REF} ${REF} 2>&1 > logs/bt2IndexRef.log
      BT2IDXCODE=$?
      if [[ $BT2IDXCODE -ne 0 ]]; then
         echo "bowtie2-build on ${REF} failed with exit code ${BT2IDXCODE}!"
         exit 4
      fi
   else
      echo "Skipping bowtie2 index creation for ${REF}"
   fi
elif [[ "${MAPPER}" == "HISAT2" ]]; then
   if [[ ! -x "$(command -v ${HISAT2BUILD})" ]]; then
      echo "hisat2-build appears to be missing, could not find at HISAT2BUILD=${HISAT2BUILD}."
      exit 13;
   fi
   #Create a hisat2 index of the reference FASTA:
   if [[ ! -e ${REF}.1.ht2 ]]; then
      echo "Making hisat2 index of ${REF}"
      ${HISAT2BUILD} -f ${REF} ${REF} 2>&1 > logs/hisat2IndexRef.log
      HISAT2IDXCODE=$?
      if [[ $HISAT2IDXCODE -ne 0 ]]; then
         echo "hisat2-build on ${REF} failed with exit code ${HISAT2IDXCODE}!"
         exit 4
      fi
   else
      echo "Skipping hisat2 index creation for ${REF}"
   fi
elif [[ "${MAPPER}" == "MM2" ]]; then
   if [[ ! -x "$(command -v ${MM2})" ]]; then
      echo "minimap2 appears to be missing, could not find at MM2=${MM2}."
      exit 12;
   fi
   #Create a minimap2 index of the reference FASTA:
   if [[ ! -e ${REF}.mmi ]]; then
      ${MM2} -d ${REF}.mmi ${REF} 2>&1 > logs/minimap2IndexRef.log
      MM2IDXCODE=$?
      if [[ $MM2IDXCODE -ne 0 ]]; then
         echo "minimap2 indexing on ${REF} failed with exit code ${MM2IDXCODE}!"
         exit 4
      fi
   else
      echo "Skipping minimap2 index creation for ${REF}"
   fi
else
   echo "Mapper ${MAPPER} not recognized, did you make a typo?"
   echo "Supported mappers: BWA, STAR, BT2 (bowtie2), MM2 (minimap2)"
   exit 2;
fi

if [[ ! -x "$(command -v ${SAMTOOLS})" ]]; then
   echo "SAMtools appears to be missing, could not find at SAMTOOLS=${SAMTOOLS}."
   exit 18;
fi
if [[ ! -e ${PICARD} ]]; then
   echo "Picard appears to be missing, could not find at PICARD=${PICARD}."
   exit 11;
fi

#Create an index of the reference FASTA:
if [[ ! -e ${REF}.fai ]]; then
   echo "Making FASTA index of ${REF}"
   ${SAMTOOLS} faidx ${REF} 2>&1 > logs/samtoolsFaidxRef.log
else
   echo "Skipping .fai creation for ${REF}"
fi

#Create a sequence dictionary of the reference FASTA (for later use by GATK):
REFDICT="${REF%.*}.dict"
if [[ ! -e ${REFDICT} ]]; then
   echo "Creating sequence dictionary for ${REF}"
   java ${JAVAMEM} -jar ${PICARD} CreateSequenceDictionary REFERENCE=${REF} OUTPUT=${REFDICT} 2>&1 > logs/picardDictRef.log
else
   echo "Skipping .dict creation for ${REF}"
fi

echo "Done creating ${MAPPER} index, .dict, and .fai for ${REF}"
