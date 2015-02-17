LUMPY
=====

A probabilistic framework for structural variant discovery

Ryan M Layer, Colby Chiang, Aaron R Quinlan, and Ira M Hall. 2014.
"LUMPY: a Probabilistic Framework for Structural Variant Discovery."
Genome Biology 15 (6): R84.
[doi:10.1186/gb-2014-15-6-r84](http://dx.doi.org/10.1186/gb-2014-15-6-r84).

## Table of Contents
1. [Quick start](#quick-start)
2. [Installation](#installation)
3. [LUMPY Express](#lumpy-express): Automated breakpoint detection for standard analyses.
4. [LUMPY (traditional)](#lumpy-traditional): Flexible and customizable breakpoint detection for advanced users.

## Quick start

Download and install
```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.
```

Run LUMPY Express
```
lumpyexpress \
    -B my.bam \
    -o output.vcf
```

## Installation

##### Requirements
- Samtools ([http://www.htslib.org/](http://www.htslib.org/))
- SAMBLASTER ([https://github.com/GregoryFaust/samblaster](https://github.com/GregoryFaust/samblaster))
- Python 2.7 ([https://www.python.org/](https://www.python.org/))
    * pysam ([https://pypi.python.org/pypi/pysam](https://pypi.python.org/pypi/pysam))
    * NumPy ([http://docs.scipy.org/doc/numpy/user/install.html](http://docs.scipy.org/doc/numpy/user/install.html))

##### Install LUMPY
```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.
```

## LUMPY Express
Automated breakpoint detection for standard analyses.

```
usage:   lumpyexpress [options]
```

**Required arguments**
```
     -B FILE  coordinate-sorted BAM file(s) (comma separated)
```

**Optional arguments**
```
     -S FILE  split reads BAM file(s) (comma separated)
     -D FILE  discordant reads BAM files(s) (comma separated)
     -o STR   output [fullBam.bam.vcf]
     -x FILE  BED file to exclude
     -P       output probability curves for each variant
     -m INT   minimum sample weight for a call [4]
     -r FLOAT trim threshold [0]
     -T DIR   temp directory [./output_prefix.XXXXXXXXXXXX]
     -k       keep temporary files
     -K FILE  path to lumpyexpress.config file
                (default: same directory as lumpyexpress)
     -v       verbose
     -h       show this message
```

#### Configuration
LUMPY Express runs several external program whose paths are specified in
[scripts/lumpyexpress.config](scripts/lumpyexpress.config). This config
must reside in the same directory as lumpyexpress, or be specified explicitly
with the -K flag.

The installation Makefile auto-generates a lumpyexpress.config file
and places it in the "bin" directory.

#### Input
LUMPY Express expects BWA-MEM aligned BAM files as input (using the -M flag to mark
shorter alignments as secondary).
It automatically parses sample, library, and read group information using the @RG
tags in the BAM header.
Each BAM file is expected to contain exactly one sample.

The minimum input is a coordinate-sorted BAM file (-B), from which LUMPY Express
extracts splitters and discordants using SAMBLASTER before running LUMPY.
Optionally, users may supply coordinate-sorted splitter (-S) and discordant (-D)
BAM files which will bypass SAMBLASTER extraction for faster analysis.

#### Output
LUMPY Express produces a VCF file according to [VCF spec 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

#### Examples

1. Run LUMPY Express on a single sample
    ```
    lumpyexpress \
        -B my.bam \
        -o output.vcf
    ```

2. Run LUMPY Express on a single sample with pre-extracted splitters and discordants
    ```
    lumpyexpress \
        -B my.bam \
        -S my.splitters.bam \
        -D my.discordants.bam \
        -o output.vcf
    ```

3. Run LUMPY Express on jointly on multiple samples
    ```
    lumpyexpress \
        -B sample1.bam,sample2.bam,sample3.bam \
        -o output.vcf
    ```

4. Run LUMPY Express jointly on multiple samples with pre-extracted splitters and discordants
    ```
    lumpyexpress \
        -B sample1.bam,sample2.bam,sample3.bam \
        -S sample1.splitters.bam,sample2.splitters..bam,sample3.splitters.bam \
        -D sample1.discordants.bam,sample2.discordants.bam,sample3.discordants.bam \
        -o output.vcf
    ```

## LUMPY (traditional)
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

    -bedpe    bedpe_file:<bedpe file>,
              id:<sample name>,
              weight:<sample weight>
```

