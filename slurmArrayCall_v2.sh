#!/bin/bash

#This script uses a metadata file and $SLURM_ARRAY_TASK_ID to
#run a pipeline job when submitted with sbatch -a

#The arguments are:
#1) job type (i.e. iADMD, IR, HC, VCFINSNP, STAR, IRRNA, MPILEUP, or MERGE)
#2) metadata file (TSV comprised of prefix, ref, read file(s))
# 2a) metadata file for MERGE is different: (merged BAM name, component BAM list)
#3) Optional override of # cores used
#4) special parameters (e.g. no_markdup, misencoded, "Alisa_v9 5")

#Note: SLURM submission parameters are specified in the sbatch
# call to this script (e.g. --mem, -N 1, --ntasks-per-node=1,
# --cpus-per-task, -t, --qos, -J)

JOBTYPE=$1
METADATA=$2
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
      else
         PREFIX="${metadatafields[0]}"
         REF="${metadatafields[1]}"
         if [[ ! -z "${metadatafields[3]}" ]]; then
            READS="${metadatafields[2]} ${metadatafields[3]}"
         else
            READS="${metadatafields[2]}"
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
elif [[ $JOBTYPE =~ "MERGE" ]]; then
   #Params: MERGED BAMLIST
   CMD="${SCRIPTDIR}/MergeBAMs.sh ${MERGED} ${BAMLIST}"
elif [[ $JOBTYPE =~ "MPILEUP" ]]; then
   #Params: PREFIX REF CORES SPECIAL
   CMD="${SCRIPTDIR}/samtoolsVariantCall.sh ${PREFIX} ${REF}${CORES} ${SPECIAL}"
else
   echo "Unintelligible job type $JOBTYPE"
   exit 3
fi

$CMD
