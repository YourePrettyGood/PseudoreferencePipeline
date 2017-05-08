#!/bin/bash
#SBATCH -J EXAMPLEiADMD
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#Tune these as appropriate, but the default number of threads used by mapping
# is 8
#SBATCH --cpus-per-task=8
#64 GB is usually a safe bet for mapping 30x of Drosophila data
#SBATCH --mem=64000

#I just use 24 hours as a timeout because that's the lowest timeout queue
# beside 1 hour, and 1 hour isn't enough for mapping in most cases.
#SBATCH -t 24:00:00
#SBATCH --qos=1day

#Adjust the following range to match the number of lines in your metadata file:
#That is, if you have 8 samples to process, set it to `-a 1-8`
#In this example, I was processing 32 individuals, so I used `-a 1-32`
#SBATCH -a 1-32

#Adjust this path to point to where you installed the PseudoreferencePipeline:
PIPELINEDIR="[Path to the Pseudoreference Pipeline]/PseudoreferencePipeline"

#Possible JOBTYPEs:
#iADMD = Align reads to a wrapped reference, sort, and mark duplicates
#STAR = Align RNAseq reads to a wrapped reference, sort, mark duplicates, and
#       split alignments across introns
#IR = Realign reads near indels using GATK IndelRealigner -- requisite before HC
#HC = Run GATK HaplotypeCaller and perform single-sample GenotypeGVCFs to
#     output one VCF record per base in the genome
#VCFINSNP = Filter the VCF, masking bases where depth or other factors are
#           inadequate, converting homozygous alt calls to the alt base, and
#           converting heterozygous calls to the degenerate IUPAC base
JOBTYPE="iADMD"

#A tab-separated-value file containing one line per sample you want to run
#The first column is a prefix to use for all intermediate and final files
#The second column is the path to the wrapped reference genome FASTA
#The third column is the path to the first read file (typ. gzipped FASTQ)
#The fourth column is the path to the second read file (if paired-end)
METADATA="Example_metadata.tsv"

#Match this to the value you used for --cpus-per-task:
NUM_THREADS=8

${PIPELINEDIR}/slurmArrayCall_v2.sh ${JOBTYPE} ${METADATA} ${NUM_THREADS}
