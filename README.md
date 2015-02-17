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
    -B sample.bam \
    -o sample
```

## Installation

#### Requirements
- Samtools [http://www.htslib.org/](http://www.htslib.org/)
- SAMBLASTER ([https://github.com/GregoryFaust/samblaster](https://github.com/GregoryFaust/samblaster))
- Python 2.7 ([https://www.python.org/](https://www.python.org/))
    * pysam ([https://pypi.python.org/pypi/pysam](https://pypi.python.org/pypi/pysam))
    * NumPy ([http://www.numpy.org/](http://www.numpy.org/))

```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.
```