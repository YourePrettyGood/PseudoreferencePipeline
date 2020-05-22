#!/bin/bash
#SBATCH -J EXAMPLEHC
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#Tune these as appropriate, but the default number of threads used by
# HaplotypeCaller is 8
#SBATCH --cpus-per-task=8
#128 GB is usually sufficient for calling variants on 30x of Drosophila data
#SBATCH --mem=128000

#24 hours is usually enough for HaplotypeCaller given enough threads
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
JOBTYPE="HC"

#A tab-separated-value file containing one line per sample you want to run
#The first column is a prefix to use for all intermediate and final files
#The second column is the path to the wrapped reference genome FASTA
#The third column is the path to the first read file (typ. gzipped FASTQ)
#The fourth column is the path to the second read file (if paired-end)
METADATA="Example_metadata.tsv"

#Match this to the value you used for --cpus-per-task:
NUM_THREADS=8

#Special options in case you have particular needs, like:
#no_HC: Skip HaplotypeCaller and go directly to GenotypeGVCFs
#no_geno: Skip GenotypeGVCFs
#mem_#g: Set the maximum memory allowed for Java during GATK calls
#no_markdup: Use a BAM that hasn't had duplicates marked
#no_IR: Use a BAM that hasn't been indel-realigned
SPECIAL=""

${PIPELINEDIR}/localArrayCall_v2.sh ${SLURM_ARRAY_TASK_ID} ${JOBTYPE} ${METADATA} ${NUM_THREADS} ${SPECIAL}
