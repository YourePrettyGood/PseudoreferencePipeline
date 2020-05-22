#All of these variables must contain the full path to the executables
#This pipeline depends on:
#BWA (can obtain via `git clone https://github.com/lh3/bwa --recursive`)
#Samtools (compatible with versions 1.0+, may have issues with older)
#BCFtools (keep in sync with samtools version)
#Tabix (keep in sync with samtools version)
#Picard (can obtain via `git clone https://github.com/broadinstitute/picard --recursive`)
#GATK (mainly used with 3.4, all steps but VCFINSNP compatible with newer)
#seqtk (can obtain via `git clone https://github.com/lh3/seqtk`)
#STAR (can obtain via `git clone https://github.com/alexdobin/STAR --recursive`)

#For PSEUDOFASTA task, also depends on:
#BEDtools (can obtain via `git clone https://github.com/arq5x/bedtools2 --recursive`)

export BWA="/usr/local/bin/bwa"
#export BT2=""
export STAR="/usr/local/bin/STAR"
export SAMTOOLS="/usr/local/bin/samtools"
export BCFTOOLS="/usr/local/bin/bcftools"
export TABIX="/usr/local/bin/tabix"
export PICARD="/home/pfreilly/Downloads/Bioinformatics/Mapping/picard/picard.jar"
export GATK="/home/pfreilly/Downloads/Bioinformatics/VariantCalling/GATK/GATK-3.4/GenomeAnalysisTK.jar"
#For GATK4:
#export GATK="/home/pfreilly/Downloads/Bioinformatics/VariantCalling/GATK/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar"

#For PSEUDOFASTA, particularly when VCFs are from GATK:
export BEDTOOLS="/usr/local/bin/bedtools"

#For heterozygosity and fixed difference rates from pseudoreferences:
export LPDS="/home/pfreilly/Downloads/Bioinformatics/General/RandomScripts/listPolyDivSites"
export NOW="/home/pfreilly/Downloads/Bioinformatics/General/RandomScripts/nonOverlappingWindows"

#If you aren't going to use the VCFINSNP step, ignore these three lines:
export SEQTK="/usr/local/bin/seqtk"
export ALISAV9="/home/pfreilly/Downloads/Bioinformatics/VariantCalling/PseudoreferencePipeline/insnp_v9_alisa.py"
export ALISAV8="/home/pfreilly/Downloads/Bioinformatics/VariantCalling/PseudoreferencePipeline/insnp_v8_alisa.py"

#For PB mapping:
export NGMLR="/home/pfreilly/.conda/envs/Snifflesenv/bin/ngmlr"
#export BLASR=""
#export MM2=""
