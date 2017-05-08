# PseudoreferencePipeline
SLURM-based pipeline for alignment, variant calling, and generation of pseudoreference FASTAs from Illumina data (DNA- or RNAseq)

## Dependencies
1. BWA (can be obtained via `git clone https://github.com/lh3/bwa --recursive`)
1. STAR (can be obtained via `git clone https://github.com/alexdobin/STAR --recursive`)
1. Samtools (mainly used with versions 1.0+, unknown compatibility below that)
1. Picard (can be obtained via `git clone https://github.com/broadinstitute/picard --recursive`)
1. GATK (mainly tested with 3.4, and IR and HC tasks should work with newer)

## Optional dependencies (if using the VCFINSNP task)
1. seqtk (can be obtained via `git clone https://github.com/lh3/seqtk`)
1. VCF to in.snp script (may be packaged here in the future)

## Usage
This pipeline semi-automates the process of alignment and variant calling for both DNAseq and RNAseq datasets.
The modular design is such that you can diagnose problems occurring at most steps by looking at that module's log.  We also take advantage of SLURM task arrays in order to make alignment and variant calling among many samples paralellized across a cluster.

**The first task in all instances of this pipeline is to index the reference genome(s) you want to use.**
**It is CRITICALLY important that your reference genome FASTA is line-wrapped for usage with GATK.**
You can use something like the `fasta_formatter` program of the 
[FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
to wrap your reference genome.

Once wrapped, then you can simply call:
`indexDictFai.sh [wrapped FASTA]`

This will index your reference genome for use with BWA, create a sequence dictionary as used by GATK, and create a FASTA index (.fai file).

Once your reference genome(s) are indexed, you can proceed with the main pipeline.

The main wrapper workhorse is slurmArrayCall_v2.sh, which gets called by your SBATCH script, and performs the specified job on a sample specified by a line in a metadata TSV file.

Take a look at the SBATCH_template_*.sbatch files for examples of how to make your SBATCH submission scripts to run this pipeline.

The metadata TSV file consists of 3 or 4 columns per line:
1. A prefix to use for all intermediate and final output files
1. The path to the **wrapped** reference genome FASTA
1. The first read file for the sample (typically a gzipped FASTQ)
1. The second read file for the sample (if it was sequenced paired-end)

Note that you can mix single-end and paired-end samples in the metadata file.

The DNAseq pipeline goes as follows:
1. Index your reference genome
1. Align using the `iADMD` task
1. Realign around indels using the `IR` task
1. Call variants using the `HC` task

The RNAseq pipeline goes as follows:
1. Create a sequence dictionary and .fai file for your reference genome
1. Align using the `STAR` task
1. Realign around indels using the `IRRNA` task
1. Call variants using the `HC` task

Note that for the RNAseq pipeline, each instance generates its own genome index for STAR, but does not generate a .dict or .fai file, so it is still necessary to perform `indexDictFai.sh` before starting the pipeline.
