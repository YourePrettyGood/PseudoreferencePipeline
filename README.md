# PseudoreferencePipeline

Pipeline for alignment, variant calling, and generation of pseudoreference FASTAs from Illumina data (DNA- or RNAseq)

The pipeline was originally written for a SLURM cluster, but has been adapted for use with GNU parallel, and can be adapted for any cluster job engine that provides an integer environment variable indicating the task ID in a task array.  For SLURM, this is $SLURM_ARRAY_TASK_ID, for SGE this is $SGE_TASK_ID, for PBS this is $PBS_ARRAYID, for LSF this is $LSB_JOBINDEX.  The pipeline has only been officially tested with SLURM and on a standalone computer, so please let me know if you have success (or problems) running it with other job engines!

## Core Dependencies

1. BWA (can be obtained via `git clone https://github.com/lh3/bwa --recursive`, must support `bwa mem`)
1. STAR (can be obtained via `git clone https://github.com/alexdobin/STAR --recursive`)
1. Samtools (mainly used with versions 1.0+, unknown compatibility below that)
1. HTSlib (a dependency of Samtools, so usually the version should match that of Samtools)
1. BCFtools (for variant calling, must be used with versions 1.7+ for `PSEUDOFASTA` task)
1. Picard (can be obtained via `git clone https://github.com/broadinstitute/picard --recursive`)
1. GATK (mainly tested with 3.4, and IR and HC tasks should work with other version 3, but not 4)
1. Unix tools (i.e. bash, GNU awk, GNU sort)

## Optional dependency (if a job engine is not an option for you)

1. GNU Parallel (if a job engine is not an option for you)

## Dependency for PSEUDOFASTA task

1. BEDtools (tested with 2.23.0+, should work with 2.3.0+)

## Installation

The scripts in this pipeline are effectively standalone, so the major installation steps are:

1. Install necessary dependencies (and remember the paths to their executables and jars)
1. Edit `pipeline_environment.sh` with the absolute paths to each of these executables and jars

In particular, the following variables should be defined in your `pipeline_environment.sh`:

For DNA-seq mapping (`MAP` task, also `DEPTH` task):
1. BWA (point to the `bwa` executable)
1. SAMTOOLS (point to the `samtools` executable)
1. PICARD (point to the `picard.jar` file)
(Will add BT2, BT2B, and MM2 for `bowtie2`, `bowtie2-build`, and `minimap2` in the future)

For RNA-seq mapping (`STAR` task):
1. STAR (point to the `STAR` executable)
(Will add HS2 and HS2B for `hisat2` and `hisat2-build` in the future)

For indel realignment and GATK variant calling (both DNA-seq and RNA-seq) (`IR`, `IRRNA`, and `HC` tasks):
1. GATK (point to the `GenomeAnalysisTK.jar` file from GATK 3.x)

For BCFtools variant calling (`MPILEUP` task):
1. BCFTOOLS (point to the `bcftools` executable)
1. TABIX (point to the `tabix` executable)

For variant filtering and site masking (`PSEUDOFASTA` task):
1. BEDTOOLS (point to the `bedtools` executable)

For calculating heterozygosity and fixed difference rates from pseudoreferences (`POLYDIV` task):
1. LPDS (point to the `listPolyDivSites` executable from [RandomScripts](YourePrettyGood/RandomScripts))
1. NOW (point to the `nonOverlappingWindows` executable from [RandomScripts](YourePrettyGood/RandomScripts))

For PacBio read mapping (`PBMAP` task):
1. NGMLR (point to the `ngmlr` executable)
(Will add MM2 for `minimap2` in the future)

## Works in progress

1. Integration of FreeBayes (depends on GNU Parallel and a bit of fancy vcflib path-work)
1. Integration of other mappers (bowtie2, minimap2, HISAT2)
1. Compatibility with GATK 4

## Usage

This pipeline semi-automates the process of alignment and variant calling for both DNAseq and RNAseq datasets.
The modular design is such that you can diagnose problems occurring at most steps by looking at that module's log.  We also take advantage of SLURM task arrays in order to make alignment and variant calling among many samples paralellized across a cluster.

Note that you need to set paths to the executables (and jar files) in pipeline_environment.sh, as this pipeline actively ignores your PATH variable.  This is intentional to ensure you know precisely which version of each dependency you are using.  Some cluster sysadmins like to change things up on you, and keep you on your toes ;)

A typical DNAseq variant calling job involves the following tasks in sequence:
1. `indexDictFai.sh` on the wrapped reference genome
1. `MAP`, aka `iADMD` (mapping, sorting, and marking of duplicates)
1. `MERGE` (optional, depending on dataset, adjusts ReadGroups and merges BAMs)
1. `IR` (indel realignment with GATK IndelRealigner)
1. `HC` or `MPILEUP` (variant calling with GATK HaplotypeCaller or BCFtools mpileup)
1. `DEPTH` (calculate windowed, per-scaffold, or genome-wide depth)
1. `PSEUDOFASTA` (filtering of variants and masking of sites to make a diploid pseudoreference FASTA)

The RNAseq pipeline involves the following tasks in sequence:
1. `indexDictFai.sh` on the wrapped reference genome
1. Align using the `STAR` task
1. (Optional if multiple libraries per sample) Merge BAMs from the `STAR` task together using the `MERGE` task
1. Realign around indels using the `IRRNA` task
1. Call variants using the `HC` task

**The first task in all instances of this pipeline is to index the reference genome(s) you want to use.**
**It is CRITICALLY important that your reference genome FASTA is line-wrapped for usage with GATK.**
You can use something like the `fasta_formatter` program of the 
[FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
to wrap your reference genome. Alternatively, there's `reformat.sh` from Brian Bushnell's 
[BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)
among other similar tools. Usage for both:

`fasta_formatter -i [original reference FASTA] -w60 -o [wrapped reference FASTA]`

or

`reformat.sh in=[original reference FASTA] -out=[wrapped reference FASTA] fastawrap=60`

Once wrapped, then you can simply call:
`[path to PseudoreferencePipeline]/indexDictFai.sh [wrapped FASTA]`

This will index your reference genome for use with BWA by default, create a sequence dictionary as used by GATK, and create a FASTA index (.fai file).

If you're not using BWA for mapping, the full usage for `indexDictFai.sh` is:

`[path to PseudoreferencePipeline]/indexDictFai.sh [wrapped FASTA] [mapper] [annotation]`

Values for `[mapper]` with tested support include `BWA` and `STAR`. Other untested values include `BT2` (bowtie2), `MM2` (minimap2), and `HISAT2` (HISAT2), however these mappers haven't been implemented in the `MAP`, `PBMAP`, and `STAR` tasks yet.

The `[annotation]` option is currently untested, but intended for future support of annotation-aware RNAseq mapping (e.g. with STAR or HISAT2).

Once your reference genome(s) are indexed, you can proceed with the main pipeline.

## Support for cluster job engines (SLURM, SGE, LSF, PBS, etc.):

The pipeline was originally developed for SLURM (hence `slurmArrayCall_v2.sh`), but generically supports any job engine that supports job arrays and an environment variable containing an integer representing the task number in the array. The tasks within the array should be numbered starting at 1, since their task ID corresponds to the line of the metadata file containing that sample.

For SLURM, the general form for any pipeline calls is:

`[path to PseudoreferencePipeline]/localArrayCall_v2.sh ${SLURM_ARRAY_TASK_ID} [TASK] [metadata file path] [number of threads] [special options]`

Be sure to specify `-a 1-[number of samples]` in your `sbatch` call for this to work. For example, for 16 lines in your metadata file, use `-a 1-16`.

Thanks to Clair Han, the pipeline has also been tested on LSF. The general form for LSF is:

`[path to PseudoreferencePipeline]/localArrayCall_v2.sh ${LSB_JOBINDEX} [TASK] [metadata file path] [number of threads] [special options]`

Be sure to specify `-J jobname[1-[number of samples]]` in your `bsub` call for this to work. For example, for 16 lines in your metadata file, use `-J myjob[1-16]`.

If anyone else would be willing to test the pipeline on their own cluster using a different job engine, please let me know (e.g. with a Github issue) so we can add support here!

## Use of the pipeline on a single machine (i.e. with GNU Parallel):

This pipeline has been adapted for use with GNU Parallel via the generic specification of task ID to `localArrayCall_v2.sh`.

An example `MAP` call using 8 cores per mapping job for a 32-line metadata file performing at most 4 jobs simultaneously might look like this:

`parallel -j4 --eta '[path to PseudoreferencePipeline]/localArrayCall_v2.sh {1} MAP example_metadata.tsv 8 2> MAP_example_job{1}.stderr > MAP_example_job{1}.stdout' ::: {1..32}`

## How to run the parallel parts of the PseudoreferencePipeline:

The main wrapper workhorse is `localArrayCall_v2.sh`, which gets called by your SBATCH script, and performs the specified job on a sample specified by a line in a metadata TSV file.

Take a look at the `SBATCH_template_*.sbatch` files for examples of how to make your SBATCH submission scripts to run this pipeline.

The typical mapping, IR, and variant calling metadata TSV file consists of 3 or 4 columns per line:
1. A prefix to use for all intermediate and final output files
1. The path to the **wrapped** reference genome FASTA
1. The first read file for the sample (typically a gzipped FASTQ)
1. The second read file for the sample (if it was sequenced paired-end)

Note that you can mix single-end and paired-end samples in the metadata file.

## Merging BAMs from multiple libraries of the same sample

If you need to merge BAMs from multiple libraries, you will need multiple metadata files. The metadata file for the `MERGE` task follows a different format, with columns:
1. Full path and name for the merged BAM file
1. Full path to the first library's BAM file (from `MAP`)
1. Full path to the second library's BAM file (from `MAP`)
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

### Generating a pseudoreference from joint genotyping:

To use `PSEUDOFASTA` with the output of `JOINTGENO` (described below), be sure to add a fourth column to your `PSEUDOFASTA` metadata file that contains the prefix of the jointly-genotyped VCF (i.e. the first column of your `JOINTGENO` metadata file.

Also, make **absolutely** sure to use the `jointgeno` SPECIAL option, otherwise `PSEUDOFASTA` will just ignore the joint prefix, and run in single-sample mode.

## Joint genotyping:

If you wish to perform joint genotyping on your samples, you'll need to use the `JOINTGENO` task. The metadata file for this task has a slightly different format than usual, with columns:
1. The prefix to use for the output jointly-genotyped VCF
1. The reference FASTA (same as the usual second column)
1. The prefix used for the previous task on the first sample
1. The prefix used for the previous task on the second sample
1. etc.

Aside from that, `JOINTGENO` just needs to know which variant caller to use as a SPECIAL option (i.e. `HC` or `MPILEUP`). 

One **major** difference between applying `JOINTGENO` to GATK HaplotypeCaller versus BCFtools is that we require that you run the `HC` task with the `no_geno` option prior to running `JOINTGENO`, whereas you can just run the BCFtools version of `JOINTGENO` from BAMs, no need to run the `MPILEUP` task beforehand. This is because joint genotyping for GATK happens at the GenotypeGVCFs step, not at the HaplotypeCaller step, so it operates on VCFs, and you will want to parallelize your HaplotypeCaller calls, which `JOINTGENO` wouldn't be able to do efficiently by itself.

Default INFO and FORMAT tags output in the jointly-genotyped VCF:
For bcftools versions < 1.3:
DP, DP4, SP

For bcftools versions >= 1.3 and < 1.7:
DP, AD, SP, ADF, ADR

For bcftools versions >= 1.7:
DP, AD, SP, ADF, ADR

## Helper tasks `DEPTH` and `POLYDIV`:

Two extra tasks are available to provide summaries of the data: `DEPTH` and `POLYDIV`

The `POLYDIV` task requires the `listPolyDivSites` and `nonOverlappingWindows` programs from [RandomScripts](https://github.com/YourePrettyGood/RandomScripts) to have paths specified in `pipeline_environment.sh` under the variables named `LPDS` and `NOW`.

The output files for `POLYDIV` have suffix `_poly_w#kb.tsv` and `_div_w#kb.tsv`, where `#` is the window size rescaled into kb.

The `DEPTH` task runs a quick `samtools flagstat` on the BAM specified by the first column of the metadata file, plus the presence or absence of the `no_markdup` and `no_IR` flags in the `SPECIAL` argument list.  It then calculates the average depth in non-overlapping windows of size specified in the `SPECIAL` argument across each scaffold.

The windowed depth values are stored in an output file with suffix `_depth_w#kb.tsv` with `#` as described for `POLYDIV`.  The flagstat results are stored in an output file with suffix `_flagstat.log`.

Two special modes exist for `DEPTH` and `POLYDIV`:
1. Genome-wide average (set SPECIAL to `w0`), in which case the output file has suffix `_depth_genomewide.tsv`, `_poly_genomewide.tsv`, or `_div_genomewide.tsv`
1. Per-scaffold averages (set SPECIAL to `w-1` or any negative number), in which case the output file has suffix `_depth_perScaf.tsv`, `_poly_perScaf.tsv`, or `_div_perScaf.tsv`

The first of these (genome-wide depth) is useful for calculating the post-markdup depth to use for a filtering criterion in the `PSEUDOFASTA` task.

## Extra/Special options:

Several of the jobtypes/tasks have special options available to cope with variations on the standard pipeline, or quirks in GATK. Multiple options may be specified by including them in a comma-separated list (no spaces allowed in the list). Special options are as follows.

`MAP` (aka `iADMD`):
1. `only_bwa`: Omits the Picard MarkDuplicates step (make sure to use `no_markdup` with any further tasks if you use this)
1. `no_markdup`: Exactly the same thing as `only_bwa`, meant for argument continuity with later tasks
1. `only_markdup`: Assumes a sorted BAM already exists with name `[PREFIX]_sorted.bam`, skips alignment with BWA-MEM, and goes straight to Picard MarkDuplicates
1. `mem_#[mg]`: Set the maximum heap size for Java (during Picard MarkDuplicates) to #[mg], e.g. 30g for 30 GB (the default) would be mem_30g
1. `interleaved`: The single FASTQ file provided is an interleaved FASTQ of paired-end reads (so use BWA's "smart pairing" mode)
1. `comments`: Include string parts after the first space in the FASTQ header in the output BAM as comments (BWA option -C)
1. `all_alignments`: Output all found alignments (extra alignments marked as secondary) (BWA option -a)

`IR`:
1. `mem_#[mg]`: Set the maximum heap size for Java (during GATK calls) to #[mg], e.g. 30g for 30 GB (the default) would be mem_30g
1. `misencoded`: FASTQ files from older Illumina sequencers have quality scores encoded as PHRED+64, and this flag converts them to PHRED+33, the standard
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)

`STAR`:
1. `no_markdup`: Skip the Picard MarkDuplicates step
1. `no_splitN`: Skip the GATK SplitNCigarReads step
1. `intronMotif`: Add the XS tag for reads aligning across canonical junctions (nominally for compatibility with Cufflinks and StringTie)

`IRRNA`:
1. `mem_#[mg]`: Set the maximum heap size for Java (during GATK calls) to #[mg], e.g. 30g for 30 GB (the default) would be mem_30g
1. `misencoded`: FASTQ files from older Illumina sequencers have quality scores encoded as PHRED+64, and this flag converts them to PHRED+33, the standard
1. `filter_mismatching_base_and_quals`: This option drops reads where the length of the sequence and quality scores differ
1. `filter_bases_not_stored`: This option drops reads of length 0 (in case your FASTQ gets mangled during adapter or quality trimming)

`HC`:
1. `no_HC`: Skips the HaplotypeCaller step and goes directly to GenotypeGVCFs
1. `no_geno`: Skips the GenotypeGVCFs step (use this prior to joint genotyping)
1. `mem_#[mg]`: Set the maximum heap size for Java (during GATK calls) to #[mg], e.g. 30g for 30 GB (the default) would be mem_30g
1. `misencoded`: Only use this option if you skipped the IR step, and have PHRED+64 encoded FASTQ files
1. `filter_mismatching_base_and_quals`: This option drops reads where the length of the sequence and quality scores differ
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)
1. `hets_#`: Sets the expected heterozygosity rate (`-hets` for GenotypeGVCFs) to the number specified (should be a decimal greater than 0.0, and less than or equal to 1.0) -- not fully implemented yet

`MPILEUP`:
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)
1. `hets_#`: Sets the expected heterozygosity rate (`-P` for `bcftools call`) to the number specified (should be a decimal greater than 0.0, and less than or equal to 1.0)

`JOINTGENO`:
1. `HC`: Use GATK HaplotypeCaller and GenotypeGVCFs for joint genotyping
1. `MPILEUP`: Use samtools/BCFtools mpileup and BCFtools call for joint genotyping
1. `no_markdup`: Uses a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)
1. `mem_#[mg]`: Set the maximum heap size for Java (during GATK calls) to #[mg], e.g. 30g for 30 GB (the default) would be mem_30g
1. `hets_#`: Sets the expected heterozygosity rate (`-hets` for `GenotypeGVCFs`, or `-P` for `bcftools call`) to the number specified (should be a decimal greater than 0.0, and less than or equal to 1.0)

`PSEUDOFASTA`:
1. `HC`: Use a VCF derived from GATK HaplotypeCaller (either the `HC` task or `JOINTGENO` with `HC`)
1. `MPILEUP`: Use a VCF derived from samtools/BCFtools mpileup and BCFtools call (either the `MPILEUP` task or `JOINTGENO` with `MPILEUP`)
1. `jointgeno`: Triggers using a jointly-genotyped VCF as input, using the first column of the metadata file as the sample ID in the VCF
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
1. `no_markdup`: Uses a pseudoref derived from a BAM in which no duplicates have been marked (e.g. `[PREFIX]_sorted.bam` or `[PREFIX]_realigned.bam`)
1. `no_IR`: Uses a pseudoref derived from a BAM that has not been indel-realigned (e.g. `[PREFIX]_sorted_markdup.bam` or `[PREFIX]_sorted.bam`)
1. `w#`: Specifies the window size (in bp) to use for calculation of non-overlapping windowed heterozygosity and fixed difference rate (usually a multiple of 1000)
1. `jointgeno`: Uses a pseudoref derived from the joint genotyping pipeline

## Typical errors

1. BAMs are missing @SQ headers and/or are corrupt
	* This is sometimes due to corrupted index files, which may be due to a race condition when indexing the reference. Please delete the index files (for the BWA pipeline, something like `rm [path to ref].fasta.* [path to ref].dict`), and run `indexDictFai.sh` **once** on the reference **prior** to starting your `iADMD` or `MAP` jobs.
	* Otherwise, run `samtools quickcheck -vvv [BAM file]` and let me know via a Github issue and/or e-mail

1. GATK RealignerTargetCreator fails
	* Most often, this is because you forgot to wrap your reference FASTA, so go back and wrap it, run `indexDictFai.sh`, and then stick the name of the wrapped FASTA in your metadata file, and it should work
	* Alternatively, this may be due to your FASTQ having PHRED+64 encoded quality scores.  If this is the case, rerun `IR` with the `misencoded` option
	* If the error has something to do with zero-length read stored in the BAM, you may have trimmed your reads in such a way as to produce 0-length reads (and not throw them out)
	* If GATK errors saying "Could not find walker with name: RealignerTargetCreator", then you need to put the path to Java 8 in your `PATH` environment variable, as you are probably using Java 11 or newer. Thanks to Clair Han for reporting and resolving this issue.

1. GATK HaplotypeCaller 
        * Occasionally, you may get an error about a read in the active region not being found in the BAM. This is typically due to multi-threading bugs in HaplotypeCaller, so try running with a different number of cores, or if you can afford it, run single-threaded. A single-threaded run should be almost guaranteed to fix this problem.
	* Work in progress, there are plenty

## Extraneous tasks

I've added a `PBMAP` task that maps PacBio data using [NGMLR](https://github.com/philres/ngmlr), but no real special options or configurations. No indexing is required, as the command line for ngmlr includes an option for not writing an index to disk. This task could be quickly expanded to mapping with [minimap2](https://github.com/lh3/minimap2), and in principle also BWA-MEM (although modern best practices seem to indicate NGMLR and minimap2 perform better). No marking of duplicates is performed, and the reads must be input as FASTQ, not as unmapped BAM.

## Future wish-list

[*] Facilitate joint genotyping (currently supported with some testing for both GATK and BCFtools)

[] Integrate FreeBayes for variant calling and `PSEUDOFASTA`
   At the moment, this is on hold because the `--report-monomorphic` option for reporting invariant sites is *extremely* slow even when parallelized, so testing is prohibitively slow.

[] Separate filtering and incorporation of indels into pseudoreferences
   This is on hold for several reasons: Many tools for downstream analysis of pseudoreferences expect them to be in the same coordinate space which incorporating indels would violate; we haven't done thorough simulations and analysis of error rates in calling indels, so filtering criteria haven't been established. It's not too hard to incorporate into the pipeline though, it just means an extra filtering/selecting step, a possible VCF merging step, and adjusting the consensus call (not sure how GATK FastaAlternateReferenceMaker handles indels though).

[] Support various other common mappers for `MAP` and `STAR`, e.g. Bowtie2, HISAT2, minimap2

[] Add compatibility with GATK 4

[] Add support for parallelizing further by genomic regions (not easy)
