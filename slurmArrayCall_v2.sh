#!/bin/bash

#This script uses a metadata file and $SLURM_ARRAY_TASK_ID to
#run a pipeline job when submitted with sbatch -a

#The arguments are:
#1) job type (i.e. iADMD, IR, HC, VCFINSNP, PSEUDOFASTA, STAR, IRRNA, MPILEUP, MERGE, POLYDIV, or DEPTH)
#2) metadata file (TSV comprised of prefix, ref, read file(s))
# 2a) metadata file for MERGE is different: (merged BAM name, component BAM list)
#3) Optional override of # cores used
#4) special parameters (e.g. no_markdup, misencoded, "Alisa_v9 5")

#Note: SLURM submission parameters are specified in the sbatch
# call to this script (e.g. --mem, -N 1, --ntasks-per-node=1,
# --cpus-per-task, -t, --qos, -J)

JOBTYPE=$1
METADATA=$2
if [[ ! -e "$METADATA" ]]; then
   echo "The metadata file does not exist! Did you make a typo? The file you specified is: ${METADATA}"
   exit 2;
fi
if [[ ! -z "$3" ]]; then
   CORES=" $3"
fi
SPECIAL=$4

WHICHSAMPLE=1
while read -r -a metadatafields
   do
   if [[ $WHICHSAMPLE -eq $SLURM_ARRAY_TASK_ID ]]; then
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
done < $2
if [[ -z "$PREFIX" && -z "$MERGED" ]]; then
   echo "Unable to find sample $SLURM_ARRAY_TASK_ID in metadata file. Skipping."
   exit 4
fi

SCRIPTDIR=`dirname $0`

if [[ $JOBTYPE =~ "iADMD" ]]; then
   #Params: PREFIX REF READS CORES
   CMD="${SCRIPTDIR}/indexAlignDictMarkDup_v2.sh ${PREFIX} ${REF} ${READS}${CORES}"
elif [[ $JOBTYPE =~ "STAR" ]]; then
   #Params: 
   CMD="${SCRIPTDIR}/STAR2passDictMarkDupSplitN_v2.sh ${PREFIX} ${REF} ${READS}${CORES}"
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
