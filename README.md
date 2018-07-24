# PseudoreferencePipeline
SLURM-based pipeline for alignment, variant calling, and generation of pseudoreference FASTAs from Illumina data (DNA- or RNAseq)

The pipeline has recently been adapted for use with GNU parallel, though this is less well-tested.

## Dependencies
1. BWA (can be obtained via `git clone https://github.com/lh3/bwa --recursive`)
1. STAR (can be obtained via `git clone https://github.com/alexdobin/STAR --recursive`)
1. Samtools (mainly used with versions 1.0+, unknown compatibility below that)
1. BCFtools (for variant calling, must be used with versions 1.7+)
1. Picard (can be obtained via `git clone https://github.com/broadinstitute/picard --recursive`)
1. GATK (mainly tested with 3.4, and IR and HC tasks should work with newer)

## Optional dependency (if SLURM is not an option for you)
1. GNU Parallel (if SLURM is not an option for you)

## Optional dependency (for PSEUDOFASTA task)
1. BEDtools (tested with 2.23.0+, should work with 2.3.0+)

## Optional dependencies (for VCFINSNP task, which is deprecated)
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

This will index your reference genome for use with BWA and STAR, create a sequence dictionary as used by GATK, and create a FASTA index (.fai file).

Once your reference genome(s) are indexed, you can proceed with the main pipeline.

The main wrapper workhorse is slurmArrayCall_v2.sh, which gets called by your SBATCH script, and performs the specified job on a sample specified by a line in a metadata TSV file.

Note: In the case when you don't have SLURM available, you can use localArrayCall_v2.sh in a call to GNU Parallel (see the end of this README for details). Hypothetically, localArrayCall_v2.sh could work for array submissions for any job engine, including SLURM, as long as an environment variable is available to indicate the current task number (like SLURM_ARRAY_TASK_ID for SLURM).  However, only GNU parallel has been tested with localArrayCall_v2.sh.

Take a look at the `SBATCH_template_*.sbatch` files for examples of how to make your SBATCH submission scripts to run this pipeline.

The metadata TSV file consists of 3 or 4 columns per line:
1. A prefix to use for all intermediate and final output files
1. The path to the **wrapped** reference genome FASTA
1. The first read file for the sample (typically a gzipped FASTQ)
1. The second read file for the sample (if it was sequenced paired-end)

Note that you can mix single-end and paired-end samples in the metadata file.

The DNAseq pipeline goes as follows:
1. Index your reference genome
1. Align using the `iADMD` task
1. (Optional if multiple libraries per sample) Merge BAMs from the `iADMD` task together using the `MERGE` task
1. Realign around indels using the `IR` task
1. Call variants using the `HC` or `MPILEUP` task

The RNAseq pipeline goes as follows:
1. Create a sequence dictionary and .fai file for your reference genome
1. Align using the `STAR` task
1. (Optional if multiple libraries per sample) Merge BAMs from the `STAR` task together using the `MERGE` task
1. Realign around indels using the `IRRNA` task
1. Call variants using the `HC` task

The final task, `PSEUDOFASTA` (replaces `VCFINSNP`, which is deprecated), filters the VCF produced by the above pipelines (though more typically the DNAseq pipeline), updates reliable variant call sites (heterozygous sites are updated to the appropriate IUPAC degenerate base), and masks unreliable sites. For the `PSEUDOFASTA` task, the fourth argument to `slurmArrayCall_v2.sh` should be a comma-separated list including the ID of the variant caller used, so either `HC` or `MPILEUP`, and special options like `no_markdup` or `no_IR` to indicate which VCF to use. The metadata TSV also has a slightly different format for the `PSEUDOFASTA` task:

Column 3 should be a filtering expression that evaluates to TRUE when a site **should** be masked. This filtering expression will be inverted to determine the sites kept/updated. The string should be a [JEXL expression](https://software.broadinstitute.org/gatk/documentation/article.php?id=1255) in the case of the SPECIAL option being `HC`, or a [BCFtools-style expression](https://samtools.github.io/bcftools/bcftools.html) in the case of `MPILEUP`.

For GATK, such a JEXL expression might look like:

`MQ <= 50.0 || DP <= 5`

In this example, sites with MQ <= 50.0 or DP <= 5 would be masked, and sites with MQ > 50.0 **and** DP > 5 would be used for updating.

A similar if not equivalent BCFtools-style expression would be:

`INFO/MQ <= 50.0 || INFO/DP <= 5`

## Merging BAMs from multiple libraries of the same sample

If you need to merge BAMs from multiple libraries, you will need multiple metadata files. The metadata file for the `MERGE` task follows a different format, with columns:
1. Full path and name for the merged BAM file
1. Full path to the first library's BAM file (from `iADMD`)
1. Full path to the second library's BAM file (from `iADMD`)
1. etc.

Be sure that the name for your merged BAM file follows the standard convention of this pipeline:
`[PREFIX]_sorted_markdup.bam`

Then you can make a new metadata file for downstream tasks (`IR`, `HC`, `MPILEUP`, etc.) in the same style as for `iADMD`, but with the first column equal to the `PREFIX` of your merged BAM.
The FASTQ columns of the metadata file are not actively used for these downstream tasks, but are checked for existence, so just substitute in the paths used for the first library that was merged.

## Extra/Special options:

Several of the jobtypes/tasks have special options available to cope with variations on the standard pipeline, or quirks in GATK. Multiple options may be specified by including them in a comma-separated list (no spaces allowed in the list). Special options are as follows.
`IR`:
1. `misencoded`: FASTQ files from older Illumina sequencers have quality scores encoded as PHRED+64, and this flag converts them to PHRED+33, the standard
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)

`IRRNA`:
1. `misencoded`: Only use this option if you skipped the IR step, and have PHRED+64 encoded FASTQ files
1. `filter_mismatching_base_and_quals`: This option drops reads where the length of the sequence and quality scores differ
1. `filter_bases_not_stored`: This option drops reads of length 0 (in case your FASTQ gets mangled during adapter or quality trimming)

`HC`:
1. `misencoded`: Only use this option if you skipped the IR step, and have PHRED+64 encoded FASTQ files
1. `filter_mismatching_base_and_quals`: This option drops reads where the length of the sequence and quality scores differ
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)

`MPILEUP`:
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)

`PSEUDOFASTA`:
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)

## Use of the pipeline without SLURM (i.e. with GNU Parallel)

This pipeline has recently been adapted for use with GNU Parallel. The only major change in usage that differs with the SLURM method is that you must pass in the ID of the TASK (i.e. the number of the line of the metadata file to use) as the first argument to localArrayCall_v2.sh.

An example `iADMD` call using 8 cores per mapping job for a 32-line metadata file performing at most 4 jobs simultaneously might look like this:

`parallel -j4 --eta '[path to PseudoreferencePipeline]/localArrayCall_v2.sh {1} iADMD example_metadata.tsv 8 2> example_job{1}_iADMD.stderr > example_job{1}_iADMD.stdout' ::: {1..32}`
