#!/bin/bash
SAMPLE=$1
#SAMPLE: Prefix used for all intermediate and output files in the pipeline
REFERENCE=$2
#REFERENCE: Path to the FASTA used as a reference for mapping
CALLER=$3
#CALLER: Short name for the variant caller used
# e.g. HC for GATK HaplotypeCaller, MPILEUP for samtools mpileup and bcftools call
SPECIAL=$4
#SPECIAL: Special options indicating input files to use, e.g. no_markdup, no_IR
FILTERSTR=$5
#FILTERSTR: Expression string to use for filtering sites -- JEXL string
# used for -select option in GATK SelectVariants for masked sites, or an
# expression string for use with bcftools filter --include
#If SPECIAL is empty, bash won't parse it as $4, instead FILTERSTR will be $4
if [[ -z "$FILTERSTR" ]]; then
   FILTERSTR="${SPECIAL}"
fi

OUTPUTDIR=""
if [[ ${SAMPLE} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${SAMPLE}`/"
   SAMPLE=`basename ${SAMPLE}`
fi

NOMARKDUP=""
REALIGNED=""
#Check that the input VCF file is there:
if [[ $SPECIAL =~ "no_markdup" ]]; then
   NOMARKDUP="_nomarkdup"
fi
if [[ ! $SPECIAL =~ "no_IR" ]]; then
   REALIGNED="_realigned"
fi
if [[ $CALLER =~ "HC" ]]; then
   VCFSUFFIX="_HC_GGVCFs.vcf"
elif [[ $CALLER =~ "MPILEUP" ]]; then
   VCFSUFFIX="mpileupcall.vcf.gz"
else
   echo "Unable to determine VCF suffix for variant caller ${CALLER}"
   exit 2
fi
INPUTVCF="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}${VCFSUFFIX}"
if [[ ! -e "${INPUTVCF}" ]]; then
   echo "Unable to find input VCF ${INPUTVCF} for variant caller ${CALLER}"
   exit 3
fi

#Load the appropriate path variables for the filtering and masking tools:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

if [[ $CALLER =~ "HC" ]]; then
   #Identify sites to mask based on the FILTERSTR:
   #Note that we only consider SNPs and invariant records, no MNPs or indels
   MASKVCF="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_sitesToMask.vcf"
   java -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${INPUTVCF} -o ${MASKVCF} -selectType SNP -selectType NO_VARIATION -select "${FILTERSTR}" 2> ${OUTPUTDIR}logs/${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_GATKSelectVariants_masking.stderr > ${OUTPUTDIR}logs/${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_GATKSelectVariants_masking.stdout
   SNPMASKCODE=$?
   if [[ $SNPMASKCODE -ne 0 ]]; then
      echo "GATK SelectVariants for masked sites on ${INPUTVCF} failed with exit code ${SNPMASKCODE}"
      exit 5
   fi
   #Extract SNPs to update by inverting the FILTERSTR:
   UPDATEVCF="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_sitesToUse.vcf"
   java -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${INPUTVCF} -o ${UPDATEVCF} -selectType SNP -selectType NO_VARIATION -select "${FILTERSTR}" -invertSelect 2> ${OUTPUTDIR}logs/${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_GATKSelectVariants_updating.stderr > ${OUTPUTDIR}logs/${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_GATKSelectVariants_updating.stdout
   SNPUPDCODE=$?
   if [[ $SNPUPDCODE -ne 0 ]]; then
      echo "GATK SelectVariants for updated sites on ${INPUTVCF} failed with exit code ${SNPUPDCODE}"
      exit 6
   fi
   #Update SNPs using FastaAlternateReferenceMaker (masking doesn't work here):
   UPDATEFASTA="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_wonky_pseudoref.fasta"
   java -jar ${GATK} -T FastaAlternateReferenceMaker -R ${REFERENCE} -V ${UPDATEVCF} -o ${UPDATEFASTA} -IUPAC ${SAMPLE} -lw 60 2> ${OUTPUTDIR}logs/${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_GATKFastaAlternateReferenceMaker.stderr > ${OUTPUTDIR}logs/${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_GATKFastaAlternateReferenceMaker.stdout
   FAUPDCODE=$?
   if [[ $FAUPDCODE -ne 0 ]]; then
      echo "GATK FastaAlternateReferenceMaker for updating sites in ${UPDATEVCF} failed with exit code ${FAUPDCODE}"
      exit 7
   fi
   #Fix the scaffold names, since GATK renames them to integers...:
   RENAMEDFASTA="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_wonky_pseudoref_renamed.fasta"
   awk 'BEGIN{FS="\t";fastarecord=1;}FNR==NR{scafs[FNR]=$1;}FNR!=NR&&/^>/{print ">"scafs[fastarecord++];}FNR!=NR&&!/^>/{print;}' ${REFERENCE}.fai ${UPDATEFASTA} > ${RENAMEDFASTA}
   FARENAMECODE=$?
   if [[ $FARENAMECODE -ne 0 ]]; then
      echo "awk script for renaming GATK FARM scaffolds failed with exit code ${FARENAMECODE}"
      exit 8
   fi
   #Finally, mask sites with bedtools maskfasta:
   PSEUDOREF="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_final_pseudoref.fasta"
   ${BEDTOOLS} maskfasta -fi ${RENAMEDFASTA} -bed ${MASKVCF} -fo ${PSEUDOREF} 2> ${OUTPUTDIR}logs/${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_bedtoolsmaskfasta.stderr > ${OUTPUTDIR}logs/${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_bedtoolsmaskfasta.stdout
   MASKFASTACODE=$?
   if [[ $MASKFASTACODE -ne 0 ]]; then
      echo "bedtools maskfasta failed with exit code ${MASKFASTACODE}"
      exit 9
   fi
elif [[ $CALLER =~ "MPILEUP" ]]; then
   UPDATEVCF="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_sitesToUse.vcf.gz"
   PSEUDOREF="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}_final_pseudoref.fasta"
   
else
   echo "Unknown variant caller ${CALLER}, making pseudoref not yet supported"
   exit 4
fi

echo "Done making pseudoreference FASTA!"
