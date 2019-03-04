# PseudoreferencePipeline
SLURM-based pipeline for alignment, variant calling, and generation of pseudoreference FASTAs from Illumina data (DNA- or RNAseq)

The pipeline has recently been adapted for use with GNU parallel, thus can be adapted for any cluster job engine that provides an integer environment variable indicating the task ID in a task array.  For SLURM, this is $SLURM_ARRAY_TASK_ID, for SGE this is $SGE_TASK_ID, for PBS this is $PBS_ARRAYID, for LSF this is $LSB_JOBINDEX.  The pipeline has only been officially tested with SLURM and on a standalone computer, so please let me know if you have success (or problems) running it with other job engines!

## Core Dependencies
1. BWA (can be obtained via `git clone https://github.com/lh3/bwa --recursive`)
1. STAR (can be obtained via `git clone https://github.com/alexdobin/STAR --recursive`)
1. Samtools (mainly used with versions 1.0+, unknown compatibility below that)
1. HTSlib (a dependency of Samtools, so usually the version should match that of Samtools)
1. BCFtools (for variant calling, must be used with versions 1.7+)
1. Picard (can be obtained via `git clone https://github.com/broadinstitute/picard --recursive`)
1. GATK (mainly tested with 3.4, and IR and HC tasks should work with newer)

## Optional dependency (if SLURM is not an option for you)
1. GNU Parallel (if SLURM is not an option for you)

## Optional dependency (for PSEUDOFASTA task)
1. BEDtools (tested with 2.23.0+, should work with 2.3.0+)

## Works in progress
1. Integration of FreeBayes (depends on GNU Parallel and a bit of fancy vcflib path-work)

## Usage
This pipeline semi-automates the process of alignment and variant calling for both DNAseq and RNAseq datasets.
The modular design is such that you can diagnose problems occurring at most steps by looking at that module's log.  We also take advantage of SLURM task arrays in order to make alignment and variant calling among many samples paralellized across a cluster.

Note that you need to set paths to the executables (and jar files) in pipeline_environment.sh, as this pipeline actively ignores your PATH variable.  This is intentional to ensure you know precisely which version of each dependency you are using.  Some cluster sysadmins like to change things up on you, and keep you on your toes ;)

A typical DNAseq variant calling job involves the following tasks in sequence:
1. indexDictFai.sh on the wrapped reference genome
1. iADMD (mapping, sorting, and marking of duplicates)
1. MERGE (optional, depending on dataset, adjusts ReadGroups and merges BAMs)
1. IR (indel realignment with GATK IndelRealigner)
1. HC or MPILEUP (variant calling with GATK HaplotypeCaller or BCFtools mpileup)
1. PSEUDOFASTA (filtering of variants and masking of sites to make a diploid pseudoreference FASTA)

The RNAseq pipeline involves the following tasks in sequence:
1. indexDictFai.sh on the wrapped reference genome
1. Align using the `STAR` task
1. (Optional if multiple libraries per sample) Merge BAMs from the `STAR` task together using the `MERGE` task
1. Realign around indels using the `IRRNA` task
1. Call variants using the `HC` task

**The first task in all instances of this pipeline is to index the reference genome(s) you want to use.**
**It is CRITICALLY important that your reference genome FASTA is line-wrapped for usage with GATK.**
You can use something like the `fasta_formatter` program of the 
[FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
to wrap your reference genome.

Once wrapped, then you can simply call:
`[path to PseudoreferencePipeline]/indexDictFai.sh [wrapped FASTA]`

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

## Indel Realignment and Variant Calling for a sample

The `IR`, `IRRNA`, `HC`, and `MPILEUP` tasks all use a metadata file with format similar to the `iADMD` or `STAR` task.  If you had to `MERGE`, just make sure the last two columns for the indel realignment/variant calling metadata file point to a real set of FASTQs, or else the job will instantly fail out.  I still need to make a simple solution to avoid this, because these steps don't actually care about the input read files.

Example calls for 8 samples using GNU parallel and SLURM:

`parallel -j8 --eta '[path to PseudoreferencePipeline]/localArrayCall_v2.sh {1} IR my_IR_metadata.tsv 1 "" 2> logs/my_IR_line{1}.stderr > logs/my_IR_line{1}.stdout' ::: {1..8}`

`[path to PseudoreferencePipeline]/localArrayCall_v2.sh ${SLURM_ARRAY_TASK_ID} IR my_IR_metadata.tsv 1 "" 2> logs/my_IR_line${SLURM_ARRAY_TASK_ID}.stderr > logs/my_IR_line${SLURM_ARRAY_TASK_ID}.stdout`

Note that the `IR` and `IRRNA` steps are not multithreaded, so specifying more than 1 core won't speed anything up.  Also, versions of BCFtools older than 1.7 do not have a multithreading option, so the same applies for `MPILEUP` in that case.

## Generating a diploid pseudoreference FASTA for a sample

The final task, `PSEUDOFASTA`, filters the VCF produced by the above pipelines (though more typically the DNAseq pipeline), updates reliable variant call sites (heterozygous sites are updated to the appropriate IUPAC degenerate base), and masks unreliable sites. For the `PSEUDOFASTA` task, the fourth argument to `slurmArrayCall_v2.sh` should be a comma-separated list including the ID of the variant caller used, so either `HC` or `MPILEUP`, and special options like `no_markdup` or `no_IR` to indicate which VCF to use. The metadata TSV also has a slightly different format for the `PSEUDOFASTA` task:

Column 3 should be a filtering expression that evaluates to TRUE when a site **should** be masked. This filtering expression will be inverted to determine the sites kept/updated. The string should be a [JEXL expression](https://software.broadinstitute.org/gatk/documentation/article.php?id=1255) in the case of the SPECIAL option being `HC`, or a [BCFtools-style expression](https://samtools.github.io/bcftools/bcftools.html) in the case of `MPILEUP`.

For GATK, such a JEXL expression might look like:

`DP <= 5 || MQ <= 50.0`

In this example, sites with DP <= 5 or MQ <= 50.0 would be masked, and sites with MQ > 50.0 **and** DP > 5 would be used for updating.

A similar if not equivalent BCFtools-style expression would be:

`INFO/MQ <= 50.0 || INFO/DP <= 5`

You can simulate in order to determine an optimized set of filtering thresholds using [VariantCallingSimulations](https://github.com/YourePrettyGood/VariantCallingSimulations).  Some preliminary results indicate that the following are reasonable thresholds for within-species mapping where pi is about 1%:

MPILEUP: DP <= about 0.5 * average post-markdup depth || MQ <= 20.0 || QUAL <= 26.0

GATK: DP <= about 0.5 * average post-markdup depth || MQ <= 50.0

## Helper tasks `DEPTH` and `POLYDIV`:

Two extra tasks are available to provide summaries of the data: `DEPTH` and `POLYDIV`

The `POLYDIV` task requires the `listPolyDivSites` and `nonOverlappingWindows` programs from [RandomScripts](https://github.com/YourePrettyGood/RandomScripts) to have paths specified in `pipeline_environment.sh` under the variables named `LPDS` and `NOW`.

The output files for `POLYDIV` have suffix `_poly_w#kb.tsv` and `_div_w#kb.tsv`, where `#` is the window size rescaled into kb.

The `DEPTH` task runs a quick `samtools flagstat` on the BAM specified by the first column of the metadata file, plus the presence or absence of the `no_markdup` and `no_IR` flags in the `SPECIAL` argument list.  It then calculates the average depth in non-overlapping windows of size specified in the `SPECIAL` argument across each scaffold.

The windowed depth values are stored in an output file with suffix `_depth_w#kb.tsv` with `#` as described for `POLYDIV`.  The flagstat results are stored in an output file with suffix `_flagstat.log`.

## Extra/Special options:

Several of the jobtypes/tasks have special options available to cope with variations on the standard pipeline, or quirks in GATK. Multiple options may be specified by including them in a comma-separated list (no spaces allowed in the list). Special options are as follows.

`iADMD`:
1. `interleaved`: The single FASTQ file provided is an interleaved FASTQ of paired-end reads (so use BWA's "smart pairing" mode)
1. `comments`: Include string parts after the first space in the FASTQ header in the output BAM as comments (BWA option -C)
1. `only_markdup`: Assumes a sorted BAM already exists with name `[PREFIX]_sorted.bam`, skips alignment with BWA-MEM, and goes straight to Picard MarkDuplicates

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
1. `no_markdup`: Uses a VCF derived from a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a VCF derived from a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)
1. `indelmaskp_#`: Specifies masking of variant calls (but not invariant sites) within # bp of either side of an indel (e.g. `indelmaskp_8` masks 8 bp on either side of an indel)
1. `indelmask_#`: Specifies masking of variant and invariant sites within # bp of either side of an indel (e.g. `indelmask_8` masks 8 bp on either side of an indel)

Note that for the `MPILEUP` caller, `indelmaskp_#` acts like adding the `-g #` flag to `bcftools filter`, and the `HC` implementation attempts to replicate this.

Also note that the `#` for `indelmask_#` must be a positive integer (so may not be 0).  The `MPILEUP` implementation would treat this properly, but the `HC` implementation would not.

`DEPTH`:
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)
1. `w#`: Specifies the window size to use for calculation of non-overlapping windowed depth (e.g. `w100000` for 100 kb non-overlapping windows)

`POLYDIV`:
1. `no_markdup`: Uses a VCF derived from a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a VCF derived from a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)
1. `w#`: Specifies the window size (in bp) to use for calculation of non-overlapping windowed heterozygosity and fixed difference rate (usually a multiple of 1000)

## Use of the pipeline without SLURM (i.e. with GNU Parallel)

This pipeline has recently been adapted for use with GNU Parallel. The only major change in usage that differs with the SLURM method is that you must pass in the ID of the TASK (i.e. the number of the line of the metadata file to use) as the first argument to localArrayCall_v2.sh.

An example `iADMD` call using 8 cores per mapping job for a 32-line metadata file performing at most 4 jobs simultaneously might look like this:

`parallel -j4 --eta '[path to PseudoreferencePipeline]/localArrayCall_v2.sh {1} iADMD example_metadata.tsv 8 2> example_job{1}_iADMD.stderr > example_job{1}_iADMD.stdout' ::: {1..32}`

## Typical errors

1. GATK RealignerTargetCreator fails
	* Most often, this is because you forgot to wrap your reference FASTA, so go back and wrap it, run `indexDictFai.sh`, and then stick the name of the wrapped FASTA in your metadata file, and it should work
	* Alternatively, this may be due to your FASTQ having PHRED+64 encoded quality scores.  If this is the case, rerun `IR` with the `misencoded` option
	* If the error has something to do with zero-length read stored in the BAM, you may have trimmed your reads in such a way as to produce 0-length reads (and not throw them out)

1. GATK HaplotypeCaller 
        * Occasionally, you may get an error about a read in the active region not being found in the BAM. This is typically due to multi-threading bugs in HaplotypeCaller, so try running with a different number of cores, or if you can afford it, run single-threaded. A single-threaded run should be almost guaranteed to fix this problem.
	* Work in progress, there are plenty

## Future wish-list

[] Integrate FreeBayes for variant calling and `PSEUDOFASTA`
   At the moment, this is on hold because the `--report-monomorphic` option for reporting invariant sites is *extremely* slow even when parallelized, so testing is prohibitively slow.

[] Facilitate joint genotyping (currently only indirectly supported for both GATK and BCFtools)
