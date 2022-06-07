*NOTE: This is LUMPY 0.2.13 with additional changes to allow lumpyexpress to function when the main file is a CRAM, not a BAM. The splitters and discordants must still be BAM files as LUMPY itself doesn't yet support CRAM as an input. This requires the hexdump command be available.*

For questions and discussion about LUMPY please visit the forum at: 

https://groups.google.com/forum/#!forum/lumpy-discuss

[![Build Status](https://travis-ci.org/arq5x/lumpy-sv.svg?branch=master)](https://travis-ci.org/arq5x/lumpy-sv) [![Support this project by running your production jobs at BatchX](https://images.batchx.io/gh-badge-logo.svg)](https://platform.batchx.io/layer-lab/tools/lumpy-sv%2Flumpyexpress "Support this project by running your production jobs at BatchX")

LUMPY
=====

A probabilistic framework for structural variant discovery.

Ryan M Layer, Colby Chiang, Aaron R Quinlan, and Ira M Hall. 2014.
"LUMPY: a Probabilistic Framework for Structural Variant Discovery."
Genome Biology 15 (6): R84.
[doi:10.1186/gb-2014-15-6-r84](http://dx.doi.org/10.1186/gb-2014-15-6-r84).

## Table of Contents
1. [Quick start](#quick-start)
2. [Installation](#installation)
3. [LUMPY Express usage](#lumpy-express-usage): Automated breakpoint detection for standard analyses.
4. [LUMPY (traditional) usage](#lumpy-traditional-usage): Flexible and customizable breakpoint detection for advanced users.
5. [Example workflows](#example-workflows)
6. [Test data](#test-data)
7. [Troubleshooting](#troubleshooting)

## Quick start

Note that [smoove](https://github.com/brentp/smoove) is the recommended way to run `lumpy` as it collects the
best-practices of `lumpy` and associated tools and will have a shorter run-time and lower false-positive rate than
`lumpyexpress` described below. 

Download and install
```
git clone --recursive https://github.com/arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.
```

Run LUMPY Express
```
lumpyexpress \
    -B my.bam \
    -S my.splitters.bam \
    -D my.discordants.bam \
    -o output.vcf
```

## Installation

##### Requirements
- LUMPY
    * g++ compiler
    * CMake
- LUMPY Express (optional)
    * Samtools (0.1.18+) ([htslib.org/](http://www.htslib.org/))
    * SAMBLASTER (0.1.19+) ([github repo](https://github.com/GregoryFaust/samblaster))
    * Python 2.7 ([python.org/](https://www.python.org/)) with pysam (0.8.3+) and NumPy (1.8.1+)
    * sambamba ([gihub repo](https://github.com/lomereiter/sambamba))
    * gawk ([GNU project](https://www.gnu.org/software/gawk/))

##### Install

Default method to install:

```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.
```

Installing with costom zlib (gzopen64 compile error):

```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
export ZLIB_PATH="/usr/lib/x86_64-linux-gnu/"; #when /usr/lib/x86_64-linux-gnu/libz.so
make
cp bin/* /usr/local/bin/.
```


## LUMPY Express usage
Automated breakpoint detection for standard analyses.

```
usage:   lumpyexpress [options]
```

**Required arguments**
```
     -B FILE  coordinate-sorted BAM file(s) (comma separated)
     -S FILE  split reads BAM file(s) (comma separated)
     -D FILE  discordant reads BAM files(s) (comma separated)

```

**Optional arguments**
```
-o STR    output [fullBam.bam.vcf]
-x FILE   BED file to exclude
-P        output probability curves for each variant
-m INT    minimum sample weight for a call [4]
-r FLOAT  trim threshold [0]
-T DIR    temp directory [./output_prefix.XXXXXXXXXXXX]
-k        keep temporary files
-K FILE   path to lumpyexpress.config file
            (default: same directory as lumpyexpress)
-v        verbose
-h        show this message
```

#### Configuration
LUMPY Express runs several external program whose paths are specified in
[scripts/lumpyexpress.config](scripts/lumpyexpress.config). This config
must reside in the same directory as lumpyexpress, or be specified explicitly
with the -K flag.

The installation Makefile auto-generates a lumpyexpress.config file
and places it in the "bin" directory.

#### Input
LUMPY Express expects BWA-MEM aligned BAM files as input.
It automatically parses sample, library, and read group information using the @RG
tags in the BAM header.
Each BAM file is expected to contain exactly one sample.

The minimum input is a coordinate-sorted BAM file (-B), from which LUMPY Express
extracts splitters and discordants using SAMBLASTER before running LUMPY.
Optionally, users may supply coordinate-sorted splitter (-S) and discordant (-D)
BAM files which will bypass SAMBLASTER extraction for faster analysis.

#### Output
LUMPY Express produces a VCF file according to [VCF spec 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

## LUMPY (traditional) usage
Flexible and customizable breakpoint detection for advanced users.

```
usage:    lumpy [options]
```

**Options**
```
-g       Genome file (defines chromosome order)
-e       Show evidence for each call
-w       File read windows size (default 1000000)
-mw      minimum weight across all samples for a call
-msw     minimum per-sample weight for a call
-tt      trim threshold
-x       exclude file bed file
-t       temp file prefix, must be to a writeable directory
-P       output probability curve for each variant
-b       output as BEDPE instead of VCF

-sr      bam_file:<file name>,
         id:<sample name>,
       	 back_distance:<distance>,
         min_mapping_threshold:<mapping quality>,
         weight:<sample weight>,
         min_clip:<minimum clip length>,
         read_group:<string>

-pe      bam_file:<file name>,
         id:<sample name>,
         histo_file:<file name>,
         mean:<value>,
         stdev:<value>,
         read_length:<length>,
         min_non_overlap:<length>,
         discordant_z:<z value>,
         back_distance:<distance>,
         min_mapping_threshold:<mapping quality>,
         weight:<sample weight>,
         read_group:<string>

-bedpe   bedpe_file:<bedpe file>,
         id:<sample name>,
         weight:<sample weight>
```

## Example workflows

#### Pre-processing
We recommend aligning data with [SpeedSeq](https://github.com/cc2qe/speedseq), which
performs BWA-MEM alignment, marks duplicates and extracts split and discordant
read-pairs.
```
speedseq align -R "@RG\tID:id\tSM:sample\tLB:lib" \
    human_g1k_v37.fasta \
    sample.1.fq \
    sample.2.fq
```

Otherwise, data may be aligned with BWA-MEM.

```
# Align the data
bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" human_g1k_v37.fasta sample.1.fq sample.2.fq \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b - \
    > sample.bam

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 sample.bam > sample.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h sample.bam \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > sample.splitters.unsorted.bam

# Sort both alignments
samtools sort sample.discordants.unsorted.bam sample.discordants
samtools sort sample.splitters.unsorted.bam sample.splitters
```

#### Running LUMPY
LUMPY has two distinct execution alternatives. LUMPY Express is a simplified wrapper for standard analyses.
LUMPY (traditional) is more customizable, for advanced users and specialized experiments.

##### LUMPY Express
- Run LUMPY Express on a single sample with pre-extracted splitters and discordants
    ```
    lumpyexpress \
        -B sample.bam \
        -S sample.splitters.bam \
        -D sample.discordants.bam \
        -o sample.vcf
    ```

- Run LUMPY Express jointly on multiple samples with pre-extracted splitters and discordants
    ```
    lumpyexpress \
        -B sample1.bam,sample2.bam,sample3.bam \
        -S sample1.splitters.bam,sample2.splitters.bam,sample3.splitters.bam \
        -D sample1.discordants.bam,sample2.discordants.bam,sample3.discordants.bam \
        -o multi_sample.vcf
    ```

- Run LUMPY Express on a tumor-normal pair
    ```
    lumpyexpress \
        -B tumor.bam,normal.bam \
        -S tumor.splitters.bam,normal.splitters.bam \
        -D tumor.discordants.bam,normal.discordants.bam \
        -o tumor_normal.vcf
    ```

##### LUMPY (traditional)
First, generate empirical insert size statistics on each library in the BAM file
```
samtools view -r readgroup1 sample.bam \
    | tail -n+100000 \
    | scripts/pairend_distro.py \
    -r 101 \
    -X 4 \
    -N 10000 \
    -o sample.lib1.histo
```
The above script (scripts/pairend_distro.py) will display mean and stdev to screen.
For these examples we will assume the mean is 500 and the stdev is 50.

- Run LUMPY with paired-end and split-reads.
    ```
    lumpy \
        -mw 4 \
        -tt 0 \
        -pe id:sample,bam_file:sample.discordants.bam,histo_file:sample.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
        -sr id:sample,bam_file:sample.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
        > sample.vcf
    ```

- Run LUMPY on a BAM file with multiple libraries.
    ```
    lumpy \
        -mw 4 \
        -tt 0 \
        -pe id:sample,read_group:rg1,bam_file:sample.discordants.bam,histo_file:sample.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
        -pe id:sample,read_group:rg2,bam_file:sample.discordants.bam,histo_file:sample.lib2.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
        -sr id:sample,bam_file:sample.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
        > sample.vcf
    ```

- Run LUMPY on multiple samples with multiple libraries.
    ```
    lumpy \
        -mw 4 \
        -tt 0 \
        -pe id:sample1,bam_file:sample1.discordants.bam,read_group:rg1,read_group:rg2,histo_file:sample1.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
        -pe id:sample1,bam_file:sample1.discordants.bam,read_group:rg3,histo_file:sample1.lib2.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
        -pe id:sample2,bam_file:sample2.discordants.bam,read_group:rg4,histo_file:sample2.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
        -sr id:sample1,bam_file:sample1.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
        -sr id:sample2,bam_file:sample2.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
        > multi_sample.vcf
    ```

- Run LUMPY with low complexity regions excluded.  
    Heng Li provides a set of low complexity regions in the supplementary information of his paper, "Toward better
    understanding of artifacts in variant calling from high-coverage samples" at
    https://doi.org/10.1093/bioinformatics/btu356. 
    ```
    unzip btu356_Supplementary_Data.zip
    unzip btu356-suppl_data.zip
    lumpy \
        -mw 4 \
        -tt 0.0 \
        -x btu356_LCR-hs37d5.bed/btu356_LCR-hs37d5.bed \
        -pe bam_file:sample.discordants.bam,histo_file:sample.pe.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:sample,min_mapping_threshold:1 \
        -sr bam_file:sample.sr.sort.bam,back_distance:10,weight:1,id:sample,min_mapping_threshold:1 \
        > sample.exclude.vcf
    ```
- Run LUMPY with regions of very high coverage excluded.  
    We can direct lumpy to ignore certain regions by using the
    exclude region option. In this example we find and then
    exclude regions that have very high coverage. First we
    use the get_coverages.py script to find the min, max, and
    mean coverages of the the sr and pe bam files, and to
    create coverage profiles for both files.

    ```
    python ../scripts/get_coverages.py \
        sample.pe.sort.bam \
	sample.sr.sort.bam
    # sample.pe.sort.bam.coverage  min:1   max:14  mean(non-zero):2.35557521272
    # sample.sr.sort.bam.coverage  min:1   max:7   mean(non-zero):1.08945936729
    ```
    
    From this output, we will choose to exclude regions that
    have more than 10x coverage.  To create the exclude file
    we will use the get_exclude_regions.py script to create
    the exclude.bed file
    ```
    python ../scripts/get_exclude_regions.py \
        10 \
	exclude.bed \
	sample.pe.sort.bam \
	sample.sr.sort.bam
    ```

    We now rerun lumpy with the exclude (-x) option
    ```
    lumpy \
        -mw 4 \
        -tt 0.0 \
        -x exclude.bed \
        -pe bam_file:sample.discordants.bam,histo_file:sample.pe.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:sample,min_mapping_threshold:1 \
        -sr bam_file:sample.sr.sort.bam,back_distance:10,weight:1,id:sample,min_mapping_threshold:1 \
        > sample.exclude.vcf
    ```

#### Post-processing
[SVTyper](https://github.com/hall-lab/svtyper) can call genotypes on LUMPY output VCF files
using a Bayesian maximum likelihood algorithm.
```
svtyper \      
    -B sample.bam \
    -S sample.splitters.bam \
    -i sample.vcf
    > sample.gt.vcf
```

## Test data
The `test/test.sh` script executes lumpy against several simulated data sets
and compares the results to the known correct result.  The sample data sets can
be found at http://layerlab.org/lumpy/data.tar.gz.  This tar ball should be
extracted into the top-level lumpy directory.  The script `test/test.sh` checks
for the the existence of this directory before running LUMPY.

## Troubleshooting
All of the bam files that lumpy processes must be position sorted. To check if your bams are sorted correctly, use the check_sorting.py script
```
python ../scripts/check_sorting.py \
    pe.pos_sorted.bam \
    sr.pos_sorted.bam \
    pe.name_sorted.bam
# pe.pos_sorted.bam
# in order
# sr.pos_sorted.bam
# in order
# pe.name_sorted.bam
# out of order:   chr10   102292476   occurred after   chr10   102292893
```
