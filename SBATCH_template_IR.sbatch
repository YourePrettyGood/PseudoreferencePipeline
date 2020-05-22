#!/bin/bash
#SBATCH -J EXAMPLEIR
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#IndelRealigner is single-threaded, so we only ever use 1 core:
#SBATCH --cpus-per-task=1
#64 GB is usually a safe bet for realigning 30x of Drosophila data
#SBATCH --mem=64000

#I just use 24 hours as a timeout because that's the lowest timeout queue
# beside 1 hour, and 1 hour isn't enough for realigning in ~1/2 of cases.
#SBATCH -t 24:00:00
#SBATCH --qos=1day

#Adjust the following range to match the number of lines in your metadata file:
#That is, if you have 8 samples to process, set it to `-a 1-8`
#In this example, I was processing 32 individuals, so I used `-a 1-32`
#SBATCH -a 1-32

#Adjust this path to point to where you installed the PseudoreferencePipeline:
PIPELINEDIR="[Path to the Pseudoreference Pipeline]/PseudoreferencePipeline"

#Possible JOBTYPEs:
#MAP = Align reads to a wrapped reference, sort, and mark duplicates
#STAR = Align RNAseq reads to a wrapped reference, sort, mark duplicates, and
#       split alignments across introns
#IR = Realign reads near indels using GATK IndelRealigner -- do this before
#     HC and MPILEUP
#IRRNA = Perform indel realignment (a la IR) on BAM derived from STAR
#HC = Run GATK HaplotypeCaller and perform single-sample GenotypeGVCFs to
#     output one VCF record per base in the genome
#MPILEUP = Run BCFtools to perform single-sample variant calling
#PSEUDOFASTA = Filter the VCF, masking bases where indicated by filtering
#              criteria, and updating the reference with genotypes at sites
#              passing these criteria
JOBTYPE="IR"

#A tab-separated-value file containing one line per sample you want to run
#The first column is a prefix to use for all intermediate and final files
#The second column is the path to the wrapped reference genome FASTA
#The third column is the path to the first read file (typ. gzipped FASTQ)
#The fourth column is the path to the second read file (if paired-end)
METADATA="Example_metadata.tsv"

#Match this to the value you used for --cpus-per-task, which should be
# 1 for IR:
NUM_THREADS=1

#Special options in case you have particular needs, like:
#mem_#g: Set the maximum memory allowed for Java during GATK calls
#misencoded: Adjust the quality scores from PHRED+64 to PHRED+33
#            since GATK doesn't like PHRED+33
SPECIAL=""

${PIPELINEDIR}/localArrayCall_v2.sh ${SLURM_ARRAY_TASK_ID} ${JOBTYPE} ${METADATA} ${NUM_THREADS} ${SPECIAL}
