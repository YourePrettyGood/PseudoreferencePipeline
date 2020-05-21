#!/bin/bash

#For BCFtools:
#TODO: Add options for ploidy specification, probably through a
# ploidy file associated with the REF, and maybe a samples TSV
# of sample IDs and sexes.
#Is ploidy specification necessary for GATK?

#TODO: Add option for specifying "expected substitution rate" (-P)

#TODO: Add -f GQ to bcftools call to add GQ annotation

#For GATK:
#TODO: Add specification of the expected heterozygosity (-hets)
#Also for -heterozygosityStandardDeviation and -indelHeterozygosity?
#Maybe also for -stand_call_conf?

#TODO: Integrate CombineGVCFs over subsets of maybe 50 samples
# to reduce GenotypeGVCFs command size and potentially avoid weird
# GATK issues (Clair noticed issues running with 100+ -V arguments
# manually through GenotypeGVCFs)
#For now, we don't bother with CombineGVCFs

#Read in the command line arguments:
#Output usage if none supplied:
if [[ $# -eq 0 ]]; then
   printf "Usage: $0 <joint VCF prefix> <reference FASTA> <special options> <list size> [list of sample prefixes]\n"
   printf "Jointly genotypes the specified set of samples.\n"
   
   exit 1
fi
#JOINTPREFIX: Prefix for jointly genotyped VCF
JOINTPREFIX=$1
#REF: Path to reference genome FASTA used for mapping
# A .fai FASTA index is expected to exist for this reference
REF=$2
#NUMPROCS: Number of threads to use for joint genotyping
# We are hard-coding this to 1 for now
#NUMPROCS=$3
NUMPROCS=1
#SPECIAL: Comma-separated list of options/flags
#'HC' or 'MPILEUP' specifies which caller to use for joint genotyping
#'no_markdup' indicates Picard MarkDuplicates was skipped in the pipeline
#'no_IR' indicates GATK IndelRealigner was skipped in the pipeline
#'hets_#' gives the heterozygosity/substitution rate to use
#         This is generally your expected theta, default is to skip this.
#
#NUMPREFIXES: Number of samples to jointly genotype
#PREFIXLIST: Space-separated list of prefixes for BAMs or GVCFs to jointly
# genotype -- length of list should match NUMPREFIXES
if [[ $4 =~ ^[0-9]+$ ]]; then
   NUMPREFIXES=$4
   PREFIXLIST="${@:5}"
else
   SPECIAL=$4
   NUMPREFIXES=$5
   PREFIXLIST="${@:6}"
fi

OUTPUTDIR=""
if [[ ${JOINTPREFIX} =~ \/ ]]; then #If the output has a path
   OUTPUTDIR="`dirname ${JOINTPREFIX}`/"
   JOINTPREFIX=`basename ${JOINTPREFIX}`
fi

#Set the paths to the executables:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

MARKDUP="_markdup"
NOMARKDUP=""
REALIGNED="_realigned"
if [[ ${SPECIAL} =~ "HC" ]]; then
   CALLER="HC"
fi
if [[ ${SPECIAL} =~ "MPILEUP" ]]; then
   CALLER="MPILEUP"
fi
if [[ ${SPECIAL} =~ "no_markdup" ]]; then
   MARKDUP=""
   NOMARKDUP="_nomarkdup"
fi
if [[ ${SPECIAL} =~ "no_IR" ]]; then
   REALIGNED=""
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
if [[ ${CALLER} == "HC" ]]; then
   echo "Running java with -Xmx${JAVAMEM} for joint genotyping of ${JOINTPREFIX}"
fi

#Check if expected heterozygosity is specified:
#Format of SPECIAL argument is hets_[01]?[.]?[0-9]*
THETA=""
THETARE='(hets)_([01]?[.]?[0-9]*)'
if [[ ${SPECIAL} =~ $THETARE ]]; then
   THETA="${BASH_REMATCH[2]}"
fi

#Check that the executables/jars are actually there:
if [[ ${CALLER} == "HC" && ! -e ${GATK} ]]; then
   echo "GATK not found, please adjust the path in the pipeline_environment script."
   echo "pipeline_environment.sh said ${GATK}"
   exit 3
fi

if [[ ${CALLER} == "MPILEUP" && ! -x "$(command -v ${BCFTOOLS})" ]]; then
   echo "BCFtools not found, please adjust the path in the pipeline_environment script."
   echo "pipeline_environment.sh said ${BCFTOOLS}"
   exit 3
fi

#REF existence check is done in *ArrayCall_v2.sh
#Check that a .fai exists for the REF in the MPILEUP case:
if [[ ${CALLER} == "MPILEUP" && ! -e ${REF}.fai ]]; then
   echo "${REF}.fai not found, but required by mpileup"
   exit 4
fi

OUTPREFIX="${JOINTPREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}"
LOGPREFIX="${OUTPUTDIR}logs/${OUTPREFIX}"
INTPREFIX="${OUTPUTDIR}${OUTPREFIX}"

#Compile either argument string (for HC) or FOFN of VCFs (for MPILEUP)
# to jointly genotype, checking for input file existence of each:
echo "Removing any existing input FOFN"
echo "rm -f ${INTPREFIX}_inputs.fofn"
rm -f ${INTPREFIX}_inputs.fofn
INPUTVCFS=""
NUMINPUTS=0
for PREFIX in ${PREFIXLIST};
   do
   if [[ ${CALLER} == "HC" ]]; then
      GZSUFFIX=""
      HCGVCF="${PREFIX}${NOMARKDUP}${REALIGNED}_HC_ERCGVCF.vcf"
      if [[ -e ${HCGVCF}.gz ]]; then
         GZSUFFIX=".gz"
      elif [[ ! -e ${HCGVCF} ]]; then
         echo "VCF for sample ${PREFIX} could not be found: ${HCGVCF} or gzipped counterpart"
         exit 5
      fi
      INPUTVCFS="${INPUTVCFS} -V ${HCGVCF}${GZSUFFIX}"
      echo "${HCGVCF}${GZSUFFIX}" >> ${INTPREFIX}_inputs.fofn
      ((NUMINPUTS++))
   elif [[ ${CALLER} == "MPILEUP" ]]; then
      if [[ -n "${REALIGNED}" ]]; then
         BAM="${PREFIX}${MARKDUP}${REALIGNED}.bam"
      else
         BAM="${PREFIX}_sorted${MARKDUP}.bam"
      fi
      if [[ -e ${BAM} ]]; then
         echo "${BAM}" >> ${INTPREFIX}_inputs.fofn
      else
         echo "BAM for sample ${PREFIX} could not be found: ${BAM}"
         exit 5
      fi
      ((NUMINPUTS++))
   else
      echo "Unable to identify variant caller ${CALLER} for cohort ${JOINTPREFIX}"
      exit 6
   fi
done
if [[ "$NUMPREFIXES" -ne "$NUMINPUTS" ]]; then
   echo "The number of samples passed to $0 (${NUMPREFIXES}) does not match the number read from the passed list (${NUMINPUTS})"
   echo "Something got mangled in the transition."
   echo "Compare the following string to the contents of your metadata file:"
   echo "${PREFIXLIST}"
   exit 7
fi

#Add in the expected heterozygosity argument if specified:
THETAOPT=""
if [[ -n "${THETA}" ]]; then
   if [[ ${CALLER} == "MPILEUP" ]]; then
      THETAOPT="-P ${THETA}"
   elif [[ ${CALLER} == "HC" ]]; then
      THETAOPT="-hets ${THETA}"
   fi
fi

if [[ ${CALLER} == "MPILEUP" ]]; then
   BCFTOOLSVERSION=`${BCFTOOLS} --version | awk '/^bcftools/{split($2, versionarr, "-"); print versionarr[1];}'`
   VERSIONARR=(${BCFTOOLSVERSION//./ })
   if [[ "${VERSIONARR[0]}" -gt "0" ]]; then
      if [[ "${VERSIONARR[0]}" -eq "1" && "${VERSIONARR[1]}" -lt "7" ]]; then
         echo "Using samtools mpileup, since BCFtools version is < 1.7"
         if [[ "${VERSIONARR[1]}" -lt "3" ]]; then
            ANNOTATIONS="DP,DP4,SP"
         else
            ANNOTATIONS="DP,AD,SP,ADF,ADR"
         fi
         MPILEUPCMD="${SAMTOOLS} mpileup -ugf ${REF} -t ${ANNOTATIONS}"
         CALLCMD="${BCFTOOLS} call -m ${THETAOPT} -Oz -o"
      else
         echo "Using bcftools mpileup, since BCFtools version is >= 1.7"
         ANNOTATIONS="FORMAT/DP,FORMAT/AD,FORMAT/SP,FORMAT/ADF,FORMAT/ADR"
         MPILEUPCMD="${BCFTOOLS} mpileup -Ou -o - --threads ${NUMPROCS} -f ${REF} -a ${ANNOTATIONS}"
         CALLCMD="${BCFTOOLS} call --threads ${NUMPROCS} -m ${THETAOPT} -Oz -o"
      fi
   else
      echo "Your version of BCFtools (${BCFTOOLSVERSION}) is too old to use, please update it"
      exit 8
   fi
fi

if [[ ${CALLER} == "HC" ]]; then
   echo "Jointly genotyping ${NUMINPUTS} samples into ${INTPREFIX}_joint.vcf.gz using GATK GenotypeGVCFs"
   echo "java -Xmx${JAVAMEM} -jar ${GATK} -T GenotypeGVCFs -nt ${NUMPROCS} -R ${REF} -allSites ${THETAOPT} -o ${INTPREFIX}_joint.vcf.gz ${INPUTVCFS} 2> ${LOGPREFIX}_GATK_GGVCFs_joint.stderr > ${LOGPREFIX}_GATK_GGVCFs_joint.stdout"
   java -Xmx${JAVAMEM} -jar ${GATK} -T GenotypeGVCFs -nt ${NUMPROCS} -R ${REF} -allSites ${THETAOPT} -o ${INTPREFIX}_joint.vcf.gz ${INPUTVCFS} 2> ${LOGPREFIX}_GATK_GGVCFs_joint.stderr > ${LOGPREFIX}_GATK_GGVCFs_joint.stdout
   GGVCFCODE=$?
   if [[ "$GGVCFCODE" -ne "0" ]]; then
      echo "Joint genotyping into ${JOINTPREFIX} failed with exit code ${GGVCFCODE}!"
      exit 9
   fi
elif [[ ${CALLER} == "MPILEUP" ]]; then
   echo "Jointly genotyping ${NUMINPUTS} samples into ${INTPREFIX}_joint.vcf.gz using BCFtools mpileup and call"
   echo "${MPILEUPCMD} -b ${INTPREFIX}_inputs.fofn 2> ${LOGPREFIX}_samtoolsMpileupJoint.stderr | ${CALLCMD} ${INTPREFIX}_joint.vcf.gz 2>&1 > ${LOGPREFIX}_bcftoolsCallJoint.log"
   ${MPILEUPCMD} -b ${INTPREFIX}_inputs.fofn 2> ${LOGPREFIX}_samtoolsMpileupJoint.stderr | ${CALLCMD} ${INTPREFIX}_joint.vcf.gz 2>&1 > ${LOGPREFIX}_bcftoolsCallJoint.log
   MPCALLCODE=$?
   if [[ "$MPCALLCODE" -ne "0" ]]; then
      echo "Samtools Mpileup or BCFtools call failed for joint genotyping for ${JOINTPREFIX} with exit code ${MPCALLCODE}!"
      exit 9
   fi
fi

echo "Joint genotyping complete for ${JOINTPREFIX} using variant caller ${CALLER}!"
