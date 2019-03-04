#!/bin/bash

#This script uses a metadata file to run a pipeline job
# when submitted as bash or GNU parallel job

#The arguments are:
#1) Task ID (line of the metadata file to use)
#2) job type (i.e. iADMD, IR, HC, VCFINSNP, PSEUDOFASTA, STAR, IRRNA, MPILEUP, MERGE, POLYDIV, or DEPTH)
#3) metadata file (TSV comprised of prefix, ref, read file(s))
# 3a) metadata file for MERGE is different: (merged BAM name, component BAM list)
#4) Optional override of # cores used
#5) special parameters (e.g. no_markdup, misencoded, "Alisa_v9 5")

TASK_ID=$1
JOBTYPE=$2
METADATA=$3
if [[ ! -e "$METADATA" ]]; then
   echo "The metadata file does not exist! Did you make a typo? The file you specified is: ${METADATA}"
   exit 2;
fi
if [[ ! -z "$4" ]]; then
   CORES=" $4"
fi
SPECIAL=$5

WHICHSAMPLE=1
while read -r -a metadatafields
   do
   if [[ $WHICHSAMPLE -eq $TASK_ID ]]; then
      if [[ $JOBTYPE =~ "MERGE" ]]; then
         MERGED="${metadatafields[0]}"
         BAMLIST="${metadatafields[@]:1}"
      elif [[ $JOBTYPE =~ "PSEUDOFASTA" ]]; then
         PREFIX="${metadatafields[0]}"
         REF="${metadatafields[1]}"
         if [[ ! -e "${REF}" ]]; then
            echo "Reference ${REF} does not exist! Did you make a typo?"
            exit 5;
         fi
         FILTERSTR="${metadatafields[@]:2}"
      else
         PREFIX="${metadatafields[0]}"
         REF="${metadatafields[1]}"
         if [[ ! -e "${REF}" ]]; then
            echo "Reference ${REF} does not exist! Did you make a typo?"
            exit 5;
         fi
         if [[ $JOBTYPE =~ "iADMD" || $JOBTYPE =~ "STAR" ]]; then
            if [[ ! -z "${metadatafields[3]}" ]]; then
               READS="${metadatafields[2]} ${metadatafields[3]}"
               if [[ ! -e "${metadatafields[3]}" ]]; then
                  echo "Read 2 file ${metadatafields[3]} does not exist! Did you make a typo?"
                  exit 7;
               fi
            else
               READS="${metadatafields[2]}"
            fi
            if [[ ! -e "${metadatafields[2]}" ]]; then
               echo "Read 1 file ${metadatafields[2]} does not exist! Did you make a typo?"
               exit 6;
            fi
         fi
      fi
   fi
   (( WHICHSAMPLE++ ))
done < $METADATA
if [[ -z "$PREFIX" && -z "$MERGED" ]]; then
   echo "Unable to find sample $TASK_ID in metadata file. Skipping."
   exit 4
fi

SCRIPTDIR=`dirname $0`

if [[ $JOBTYPE =~ "iADMD" ]]; then
   #Params: PREFIX REF READS CORES SPECIAL
   CMD="${SCRIPTDIR}/indexAlignDictMarkDup_v2.sh ${PREFIX} ${REF} ${READS}${CORES} ${SPECIAL}"
elif [[ $JOBTYPE =~ "STAR" ]]; then
   #Params: PREFIX REF READS CORES SPECIAL
   CMD="${SCRIPTDIR}/STAR2passDictMarkDupSplitN_v2.sh ${PREFIX} ${REF} ${READS}${CORES} ${SPECIAL}"
elif [[ $JOBTYPE =~ "IRRNA" ]]; then
   #Params: PREFIX REF SPECIAL
   CMD="${SCRIPTDIR}/IndelRealignmentRNAseq.sh ${PREFIX} ${REF} ${SPECIAL}"
elif [[ $JOBTYPE =~ "IR" ]]; then
   #Params: PREFIX REF SPECIAL
   CMD="${SCRIPTDIR}/IndelRealignment.sh ${PREFIX} ${REF} ${SPECIAL}"
elif [[ $JOBTYPE =~ "HC" ]]; then
   #Params: PREFIX REF CORES SPECIAL
   if [[ ! -z "$SPECIAL" ]]; then
      if [[ -z "$CORES" ]]; then
         CORES=" 8"
      fi
   fi
   CMD="${SCRIPTDIR}/HaplotypeCaller.sh ${PREFIX} ${REF}${CORES} ${SPECIAL}"
elif [[ $JOBTYPE =~ "VCFINSNP" ]]; then
   #Params: PREFIX REF SCRIPT MINDEPTH
   CMD="${SCRIPTDIR}/vcf_to_insnp.sh ${PREFIX} ${REF} ${SPECIAL}"
elif [[ $JOBTYPE =~ "PSEUDOFASTA" ]]; then
   #Params: PREFIX REF CALLER SPECIAL FILTERSTR
   #Lazy way would be to extract CALLER from SPECIAL and that's it
   IFS="," read -r -a specialops <<< "${SPECIAL}"
   CALLER="${specialops[0]}"
   CMD="${SCRIPTDIR}/vcfToPseudoref.sh ${PREFIX} ${REF} ${CALLER} ${SPECIAL} ${FILTERSTR}"
elif [[ $JOBTYPE =~ "MERGE" ]]; then
   #Params: MERGED BAMLIST
   CMD="${SCRIPTDIR}/MergeBAMs.sh ${MERGED} ${BAMLIST}"
elif [[ $JOBTYPE =~ "MPILEUP" ]]; then
   #Params: PREFIX REF CORES SPECIAL
   CMD="${SCRIPTDIR}/samtoolsVariantCall.sh ${PREFIX} ${REF}${CORES} ${SPECIAL}"
elif [[ $JOBTYPE =~ "POLYDIV" ]]; then
   #Params: PREFIX REF SPECIAL
   CMD="${SCRIPTDIR}/polydiv.sh ${PREFIX} ${REF} ${SPECIAL}"
elif [[ $JOBTYPE =~ "DEPTH" ]]; then
   #Params: PREFIX SPECIAL
   CMD="${SCRIPTDIR}/flagstatdepth.sh ${PREFIX} ${SPECIAL}"
else
   echo "Unintelligible job type $JOBTYPE"
   exit 3
fi

$CMD
