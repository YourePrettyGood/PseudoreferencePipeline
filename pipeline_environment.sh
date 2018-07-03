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

export BWA="[Path to BWA]/bwa"
export STAR="[Path to STAR]/STAR"
export SAMTOOLS="[Path to Samtools]/samtools"
export BCFTOOLS="[Path to BCFtools]/bcftools"
export TABIX="[Path to Tabix]/tabix"
export PICARD="[Path to Picard]/picard.jar"
export GATK="[Path to GATK]/GenomeAnalysisTK.jar"

#For PSEUDOFASTA, particularly when VCFs are from GATK:
export BEDTOOLS="[Path to BEDtools]/bedtools"

#For heterozygosity and fixed difference rates from pseudoreferences:
export LPDS="[Path to listPolyDivSites]/listPolyDivSites"
export NOW="[Path to nonOverlappingWindows]/nonOverlappingWindows"

#If you aren't going to use the VCFINSNP step, ignore these three lines:
export SEQTK="[Path to seqtk]/seqtk"
export ALISAV9="[Path to PseudoreferencePipeline]/insnp_v9_alisa.py"
export ALISAV8="[Path to PseudoreferencePipeline]/insnp_v8_alisa.py"
