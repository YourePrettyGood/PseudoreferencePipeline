#!/bin/bash
SAMPLE=$1
REFERENCE=$2
NPROCS=$3
if [[ -z "$3" ]]; then
   NPROCS=8
fi
if [[ -e ${SAMPLE}_HC_ERCGVCF.vcf ]]; then
   java -Xmx30g -jar $GATK -T GenotypeGVCFs -nt ${NPROCS} -R ${REFERENCE} -V ${SAMPLE}_HC_ERCGVCF.vcf -allSites -o ${SAMPLE}_HC_GGVCFs.vcf 2>&1 > ${SAMPLE}_GATK_GGVCFs.log
else
   echo "HaplotypeCaller failed"
fi
