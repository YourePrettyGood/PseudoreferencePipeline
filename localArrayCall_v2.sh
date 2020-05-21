#!/bin/bash

#This script uses a metadata file to run a pipeline job
# when submitted as bash or GNU parallel job

#The arguments are:
#1) Task ID (line of the metadata file to use)
#2) job type (i.e. iADMD aka MAP, IR, HC, VCFINSNP, PSEUDOFASTA, STAR, IRRNA, MPILEUP, MERGE, JOINTGENO, POLYDIV, DEPTH, or PBMAP)
#3) metadata file (TSV comprised of prefix, ref, read file(s))
# 3a) metadata file for MERGE is different: (merged BAM name, component BAM list)
# 3b) metadata file for PSEUDOFASTA is different: (prefix, ref, VCF exclusive filter criteria)
# 3c) metadata file for JOINTGENO is different: (joint VCF prefix, list of component sample prefixes)
#4) Optional override of # cores used
#5) special parameters (e.g. no_markdup, misencoded, "Alisa_v9 5", etc.)

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
while IFS=$'\a' read -r -a metadatafields
   do
   if [[ $WHICHSAMPLE -eq $TASK_ID ]]; then
      if [[ $JOBTYPE =~ "MERGE" ]]; then
         if [[ ${#metadatafields[@]} -lt "2" ]]; then
            echo "No list of BAMs to merge supplied! Is your metadata file actually tab-separated (not space-separated)?"
            exit 8;
         fi
         MERGED="${metadatafields[0]}"
         BAMLIST="${metadatafields[@]:1}"
      elif [[ $JOBTYPE =~ "JOINTGENO" ]]; then
         if [[ ${#metadatafields[@]} -lt "3" ]]; then
            echo "No list of prefixes to jointly genotype supplied! Is your metadata file actually tab-separated (not space-separated)?"
            exit 8;
         fi
         JOINTPREFIX="${metadatafields[0]}"
         REF="${metadatafields[1]}"
         if [[ ! -e "${REF}" ]]; then
            echo "Reference ${REF} does not exist! Did you make a typo?"
            exit 5;
         fi
         PREFIXLIST="${metadatafields[@]:2}"
         NUMPREFIXES=${#metadatafields[@]}
         ((NUMPREFIXES-=2))
      elif [[ $JOBTYPE =~ "PSEUDOFASTA" ]]; then
         if [[ ${#metadatafields[@]} -lt "3" ]]; then
            echo "No string of VCF site filters supplied! Is your metadata file actually tab-separated (not space-separated)?"
            exit 8;
         fi
         PREFIX="${metadatafields[0]}"
         REF="${metadatafields[1]}"
         if [[ ! -e "${REF}" ]]; then
            echo "Reference ${REF} does not exist! Did you make a typo?"
            exit 5;
         fi
#         FILTERSTR="${metadatafields[@]:2}"
         FILTERSTR="${metadatafields[2]}"
         JOINTPREFIX="${metadatafields[3]}"
      else
         if [[ ${#metadatafields[@]} -lt "3" ]]; then
            echo "No read files supplied! Is your metadata file actually tab-separated (not space-separated)?"
            exit 8;
         fi
         PREFIX="${metadatafields[0]}"
         REF="${metadatafields[1]}"
         if [[ ! -e "${REF}" ]]; then
            echo "Reference ${REF} does not exist! Did you make a typo?"
            exit 5;
         fi
         if [[ $JOBTYPE =~ "iADMD" || $JOBTYPE =~ "MAP" || $JOBTYPE =~ "STAR" || $JOBTYPE =~ "PBMAP" ]]; then
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
done < <(tr "\t" "\a" < $METADATA)
if [[ -z "$PREFIX" && -z "$MERGED" && -z "${JOINTPREFIX}" ]]; then
   echo "Unable to find sample $TASK_ID in metadata file. Skipping."
   exit 4
fi

SCRIPTDIR=`dirname $0`

if [[ $JOBTYPE =~ "PBMAP" ]]; then
   #PBMAP has to be before MAP because we use regex match, and
   # PBMAP does regex match to "MAP", even though we added PBMAP
   # later in development.
   #Params: PREFIX REF READS CORES SPECIAL
   CMD="${SCRIPTDIR}/PBalign.sh ${PREFIX} ${REF} ${READS}${CORES} ${SPECIAL}"
elif [[ $JOBTYPE =~ "iADMD" || $JOBTYPE =~ "MAP" ]]; then
   #Params: PREFIX REF READS CORES SPECIAL
   CMD="${SCRIPTDIR}/indexAlignDictMarkDup_v2.sh ${PREFIX} ${REF} ${READS}${CORES} ${SPECIAL}"
elif [[ $JOBTYPE =~ "STAR" ]]; then
   #Params: PREFIX REF READS CORES SPECIAL
   CMD="${SCRIPTDIR}/STAR2passDictMarkDupSplitN_v2.sh ${PREFIX} ${REF} ${READS}${CORES} ${SPECIAL}"
elif [[ $JOBTYPE =~ "MERGE" ]]; then
   #Params: MERGED BAMLIST
   CMD="${SCRIPTDIR}/MergeBAMs.sh ${MERGED} ${BAMLIST}"
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
elif [[ $JOBTYPE =~ "MPILEUP" ]]; then
   #Params: PREFIX REF CORES SPECIAL
   CMD="${SCRIPTDIR}/samtoolsVariantCall.sh ${PREFIX} ${REF}${CORES} ${SPECIAL}"
elif [[ $JOBTYPE =~ "JOINTGENO" ]]; then
   #Params: JOINTPREFIX REF CORES SPECIAL NUMPREFIXES PREFIXLIST
   CMD="${SCRIPTDIR}/JointGenotyping.sh ${JOINTPREFIX} ${REF} ${CORES} ${SPECIAL} ${NUMPREFIXES} ${PREFIXLIST}"
elif [[ $JOBTYPE =~ "PSEUDOFASTA" ]]; then
   #Params: PREFIX REF CALLER SPECIAL JOINTPREFIX FILTERSTR
   #Lazy way would be to extract CALLER from SPECIAL and that's it
   IFS="," read -r -a specialops <<< "${SPECIAL}"
   CALLER="${specialops[0]}"
   CMD="${SCRIPTDIR}/vcfToPseudoref.sh ${PREFIX} ${REF} ${CALLER} ${SPECIAL} ${JOINTPREFIX} ${FILTERSTR}"
elif [[ $JOBTYPE =~ "POLYDIV" ]]; then
   #Params: PREFIX REF SPECIAL
   CMD="${SCRIPTDIR}/polydiv.sh ${PREFIX} ${REF} ${SPECIAL}"
elif [[ $JOBTYPE =~ "DEPTH" ]]; then
   #Params: PREFIX SPECIAL
   CMD="${SCRIPTDIR}/flagstatdepth.sh ${PREFIX} ${SPECIAL}"
elif [[ $JOBTYPE =~ "VCFINSNP" ]]; then
   #Params: PREFIX REF SCRIPT MINDEPTH
   echo "Warning: VCFINSNP task is deprecated"
   CMD="${SCRIPTDIR}/vcf_to_insnp.sh ${PREFIX} ${REF} ${SPECIAL}"
else
   echo "Unintelligible job type $JOBTYPE"
   exit 3
fi

$CMD
