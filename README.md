LUMPY
=====

A probabilistic framework for structural variant discovery

Ryan M Layer, Colby Chiang, Aaron R Quinlan, and Ira M Hall. 2014. "LUMPY: a Probabilistic Framework for Structural Variant Discovery." Genome Biology 15 (6): R84. [doi:10.1186/gb-2014-15-6-r84](http://dx.doi.org/10.1186/gb-2014-15-6-r84).

<!---
## Table of Contents
1. [Requirements](#requirements)
2. [Quick start](#quick-start)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Example workflows](#example-workflows)
-->

## Quick start
```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.

lumpyexpress \
    -B my.bam \
    -o output_prefix
```

## Installation

#### Requirements
- Samtools ([http://www.htslib.org/](http://www.htslib.org/))
- SAMBLASTER ([https://github.com/GregoryFaust/samblaster](https://github.com/GregoryFaust/samblaster))
- Python 2.7 ([https://www.python.org/](https://www.python.org/))
    * pysam ([https://pypi.python.org/pypi/pysam](https://pypi.python.org/pypi/pysam))
    * NumPy ([http://docs.scipy.org/doc/numpy/user/install.html](http://docs.scipy.org/doc/numpy/user/install.html))

```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.
```

## Usage

[LUMPY Express](#lumpy-express): automated script for standard analyses  
[LUMPY (traditional)](#lumpy-traditional): flexible and customizable for specialty analyses

#### LUMPY Express
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
     -o STR   output prefix [fullBam.bam]
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

#### LUMPY (traditional)
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




