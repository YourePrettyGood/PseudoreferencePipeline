#!/bin/bash
SAMPLE="$1"
#SAMPLE: Prefix used for all intermediate and output files in the pipeline
REFERENCE="$2"
#REFERENCE: Path to the FASTA used as a reference for mapping
CALLER="$3"
#CALLER: Short name for the variant caller used
# e.g. HC for GATK HaplotypeCaller, MPILEUP for samtools mpileup and bcftools call
SPECIAL="$4"
#SPECIAL: Special options indicating input files to use, e.g. no_markdup, no_IR
FILTERSTR="${@:5}"
echo "${FILTERSTR}"
#FILTERSTR: Expression string to use for filtering sites -- JEXL string
# used for -select option in GATK SelectVariants for masked sites, or an
# expression string for use with bcftools filter --include
#If SPECIAL is empty, bash won't parse it as $4, instead FILTERSTR will be $4
if [[ -z "$FILTERSTR" ]]; then
   FILTERSTR="${@:4}"
fi

#Check if we want to mask positives around indels:
INDELWINDOW=""
INDELMASKRE='indelmask_([0-9]+)'
if [[ ${SPECIAL} =~ $INDELMASKRE ]]; then
   INDELWINDOW="${BASH_REMATCH[1]}"
   echo "Masking ${INDELWINDOW} bp around each indel for sample ${SAMPLE}"
fi

OUTPUTDIR=""
if [[ ${SAMPLE} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${SAMPLE}`/"
   SAMPLE=`basename ${SAMPLE}`
fi

mkdir -p ${OUTPUTDIR}logs

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
   VCFSUFFIX="_mpileupcall.vcf.gz"
else
   echo "Unable to determine VCF suffix for variant caller ${CALLER}"
   exit 2
fi

OUTPREFIX="${SAMPLE}${NOMARKDUP}${REALIGNED}_${CALLER}"

INPUTVCF="${OUTPUTDIR}${SAMPLE}${NOMARKDUP}${REALIGNED}${VCFSUFFIX}"
if [[ ! -e "${INPUTVCF}" ]]; then
   if [[ $CALLER =~ "HC" && -e "${INPUTVCF}.gz" ]]; then
      INPUTVCF="${INPUTVCF}.gz"
   else
      echo "Unable to find input VCF ${INPUTVCF} for variant caller ${CALLER}"
      exit 3
   fi
fi

#Load the appropriate path variables for the filtering and masking tools:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

LOGPREFIX="${OUTPUTDIR}logs/${OUTPREFIX}"
INTPREFIX="${OUTPUTDIR}${OUTPREFIX}"

if [[ $CALLER =~ "HC" ]]; then
   #If indelmask has been specified, identify positives near indels
   INDELMASK=""
   if [[ ! -z "${INDELWINDOW}" ]]; then
      #Identify indel positions, and output their windowed regions:
      echo "Identifying indel positions and outputting windowed flanking regions for sample ${SAMPLE}"
      INDELMASKINTERVALS="${INTPREFIX}_indelmaskintervals.bed"
      java -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${INPUTVCF} -selectType INDEL 2> ${LOGPREFIX}_GATKSelectVariants_indels.stderr | ${SCRIPTDIR}/GATK_indel_windows.awk - ${INDELWINDOW} | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - 2> ${LOGPREFIX}_bedtoolsMergeIndelMask.stderr > ${INDELMASKINTERVALS}
      INDELINTCODE=$?
      if [[ $INDELINTCODE -ne 0 ]]; then
         echo "GATK SelectVariants for indels failed for sample ${SAMPLE} with exit code ${INDELINTCODE}"
         exit 14
      fi
      #Now extract positives within these intervals:
      echo "Extracting variant calls within flanking regions of indels for sample ${SAMPLE}"
      INDELMASKVCF="${INTPREFIX}_indelmask_${INDELWINDOW}.vcf"
      java -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${INPUTVCF} -o ${INDELMASKVCF} --intervals ${INDELMASKINTERVALS} -selectType SNP 2> ${LOGPREFIX}_GATKSelectVariants_indelmask_${INDELWINDOW}.stderr > ${LOGPREFIX}_GATKSelectVariants_indelmask_${INDELWINDOW}.stdout
      INDELMASKCODE=$?
      if [[ $INDELMASKCODE -ne 0 ]]; then
         echo "GATK SelectVariants for SNPs near indels failed for sample ${SAMPLE} with exit code ${INDELMASKCODE}"
         exit 15
      fi
      INDELMASK="--excludeIntervals ${INDELMASKVCF}"
   fi
   #Extract SNPs to update by inverting the FILTERSTR:
   #Be sure to exclude sites with uncalled genotypes
   #They lack INFO/DP, and GATK JEXL doesn't handle FORMAT/DP (correctly?)
   echo "Extracting SNPs and invariant sites passing the filters for sample ${SAMPLE}"
   UPDATEVCF="${INTPREFIX}_sitesToUse.vcf"
   java -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${INPUTVCF} -o ${UPDATEVCF} ${INDELMASK} -selectType SNP -selectType NO_VARIATION -select "${FILTERSTR}" -invertSelect 2> ${LOGPREFIX}_GATKSelectVariants_updating.stderr > ${LOGPREFIX}_GATKSelectVariants_updating.stdout
   SNPUPDCODE=$?
   if [[ $SNPUPDCODE -ne 0 ]]; then
      echo "GATK SelectVariants for updated sites on ${INPUTVCF} failed with exit code ${SNPUPDCODE}"
      exit 5
   fi
   #Identify sites to mask by complementing the sites to update:
   echo "Constructing BED of sites to mask for sample ${SAMPLE}"
   MASKINGBED="${INTPREFIX}_sitesToMask.bed"
   awk 'BEGIN{FS="\t";OFS="\t";}!/^#/{split($9, formatarr, ":"); for (elem in formatarr) {if (formatarr[elem] == "GT") {gtindex=elem;break;};}; split($10, samplearr, ":"); if (samplearr[gtindex] != "./.") {print $1, $2-1, $2;};}' ${UPDATEVCF} | sort -k1,1 -k2,2n -k3,3n | ${BEDTOOLS} merge -i - | ${BEDTOOLS} complement -i - -g <(cut -f1,2 ${REFERENCE}.fai) 2> ${LOGPREFIX}_bedtoolscomplement.stderr > ${MASKINGBED}
   MASKBEDCODE=$?
   if [[ $MASKBEDCODE -ne 0 ]]; then
      echo "bedtools merge or bedtools complement for sample ${SAMPLE} failed with exit code ${MASKBEDCODE}"
      exit 6
   fi
   #Update SNPs using FastaAlternateReferenceMaker (masking doesn't work here):
   echo "Updating FASTA with SNPs passing filters for sample ${SAMPLE}"
   UPDATEFASTA="${INTPREFIX}_wonky_pseudoref.fasta"
   java -jar ${GATK} -T FastaAlternateReferenceMaker -R ${REFERENCE} -V ${UPDATEVCF} -o ${UPDATEFASTA} -IUPAC ${SAMPLE} -lw 60 2> ${LOGPREFIX}_GATKFastaAlternateReferenceMaker.stderr > ${LOGPREFIX}_GATKFastaAlternateReferenceMaker.stdout
   FAUPDCODE=$?
   if [[ $FAUPDCODE -ne 0 ]]; then
      echo "GATK FastaAlternateReferenceMaker for updating sites in ${UPDATEVCF} failed with exit code ${FAUPDCODE}"
      exit 7
   fi
   #Fix the scaffold names, since GATK renames them to integers...:
   echo "Fixing FASTA scaffold names for sample ${SAMPLE} because GATK renames them to integers..."
   RENAMEDFASTA="${INTPREFIX}_wonky_pseudoref_renamed.fasta"
   awk 'BEGIN{FS="\t";fastarecord=1;}FNR==NR{scafs[FNR]=$1;}FNR!=NR&&/^>/{print ">"scafs[fastarecord++];}FNR!=NR&&!/^>/{print;}' ${REFERENCE}.fai ${UPDATEFASTA} > ${RENAMEDFASTA}
   FARENAMECODE=$?
   if [[ $FARENAMECODE -ne 0 ]]; then
      echo "awk script for renaming GATK FARM scaffolds for sample ${SAMPLE} failed with exit code ${FARENAMECODE}"
      exit 8
   fi
   #Finally, mask sites with bedtools maskfasta:
   echo "Masking sites failing the filters for sample ${SAMPLE}"
   PSEUDOREF="${INTPREFIX}_final_pseudoref.fasta"
   ${BEDTOOLS} maskfasta -fi ${RENAMEDFASTA} -bed ${MASKINGBED} -fo ${PSEUDOREF} 2> ${LOGPREFIX}_bedtoolsmaskfasta.stderr > ${LOGPREFIX}_bedtoolsmaskfasta.stdout
   MASKFASTACODE=$?
   if [[ $MASKFASTACODE -ne 0 ]]; then
      echo "bedtools maskfasta for sample ${SAMPLE} failed with exit code ${MASKFASTACODE}"
      exit 9
   fi
elif [[ $CALLER =~ "MPILEUP" ]]; then
   UPDATEVCF="${INTPREFIX}_filtered.vcf.gz"
   SITESTOUSEVCF="${INTPREFIX}_sitesToUse.vcf.gz"
   MASKINGBED="${INTPREFIX}_sitesToMask.bed"
   PSEUDOREF="${INTPREFIX}_final_pseudoref.fasta"
   #Mask positives within INDELWINDOW around each indel, if INDELWINDOW is set:
   INDELMASK=""
   if [[ ! -z "${INDELWINDOW}" ]]; then
      INDELMASK="--SnpGap ${INDELWINDOW}"
   fi
   #Filter sites using the filtering expression:
   echo "Applying filters to the input VCF for sample ${SAMPLE}"
   ${BCFTOOLS} filter -mx -S. -Oz -sFAIL -e "${FILTERSTR}" ${INDELMASK} -o ${UPDATEVCF} ${INPUTVCF} 2> ${LOGPREFIX}_bcftoolsfilter.stderr > ${LOGPREFIX}_bcftoolsfilter.stdout
   FILTERCODE=$?
   if [[ $FILTERCODE -ne 0 ]]; then
      echo "bcftools filter for sample ${SAMPLE} failed with exit code ${FILTERCODE}"
      exit 10
   fi
   #Generate BED of masked sites:
   echo "Constructing BED of filtered variant and invariant sites for sample ${SAMPLE}"
   ${BCFTOOLS} query -i 'FILTER=="PASS" && (TYPE=="SNP" || TYPE=="REF")' -f '%CHROM\t%POS0\t%POS\n' ${UPDATEVCF} 2> ${LOGPREFIX}_bcftoolsqueryused.stderr | ${BEDTOOLS} complement -i - -g <(cut -f1,2 ${REFERENCE}.fai) 2> ${LOGPREFIX}_bedtoolscomplement.stderr > ${MASKINGBED}
   MASKBEDCODE=$?
   if [[ $MASKBEDCODE -ne 0 ]]; then
      echo "bcftools query or bedtools complement for sample ${SAMPLE} failed with exit code ${MASKBEDCODE}"
      exit 11
   fi
   #Extract sites to use:
   #We skip REF sites, since they don't need updating, and that makes the VCF smaller
   echo "Extracting variant sites that pass the filters for updating for sample ${SAMPLE}"
   ${BCFTOOLS} view -i 'FILTER=="PASS" && TYPE=="SNP"' -Oz -o ${SITESTOUSEVCF} ${UPDATEVCF} 2> ${LOGPREFIX}_bcftoolsviewsitestouse.stderr > ${LOGPREFIX}_bcftoolsviewsitestouse.stdout
   SITESTOUSECODE=$?
   if [[ $SITESTOUSECODE -ne 0 ]]; then
      echo "bcftools view to extract used sites for sample ${SAMPLE} failed with error code ${SITESTOUSECODE}"
      exit 12
   fi
   echo "Indexing sitesToUse.vcf.gz for sample ${SAMPLE}"
   ${TABIX} ${SITESTOUSEVCF}
   #Generate pseudoref:
   echo "Generating pseudoreference FASTA for sample ${SAMPLE}"
   ${BCFTOOLS} consensus --iupac-codes -f ${REFERENCE} -m ${MASKINGBED} -s ${SAMPLE} -o ${PSEUDOREF} ${SITESTOUSEVCF} 2> ${LOGPREFIX}_bcftoolsconsensus.stderr > ${LOGPREFIX}_bcftoolsconsensus.stdout
   PSEUDOREFCODE=$?
   if [[ $PSEUDOREFCODE -ne 0 ]]; then
      echo "bcftools consensus for sample ${SAMPLE} failed with exit code ${PSEUDOREFCODE}"
      exit 13
   fi
else
   echo "Unknown variant caller ${CALLER}, making pseudoref not yet supported"
   exit 4
fi

echo "Done making pseudoreference FASTA for sample ${SAMPLE} and variant caller ${CALLER}!"
