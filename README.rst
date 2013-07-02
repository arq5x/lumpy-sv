Overview
========

This repository contains the source code for the Lumpy structural variation
detection framework developed in the Hall and Quinlan laboratories at the
University of Virginia.

Installation
============
Lumpy installation requires the GNU Scientific Library (GSL). We also recommend
installing samtools, bedtools, bamtools, novoalign (or bwa), yaha, and
Statistics::Descriptive.  Below is a step-by-step tutorial for how to install
and use Lumpy. This guide assumes that `/usr/local/bin` is writable and in your
path.  If either is not true, use another directory that is both writable and
in your path, or contact your administrator.  If you have questions, email me.

Install the GNU Scientific Libraries (GSL).

    - This is typically quite simple, as one can use package managers.
        - e.g., for OS X using Homebrew: `brew install gsl`
        - e.g., for Ubuntu: `apt-get install gsl`

Install samtools
::
    git clone git://github.com/samtools/samtools.git
    cd samtools
    make
    cp samtools /usr/local/bin/.

Install bedtools
::
    git clone git@github.com:arq5x/bedtools.git    
    cd bedtools
    make
    cp bin/bedtools /usr/local/bin/.

Install bamtools
::
    git clone git://github.com/pezmaster31/bamtools.git
    cd bamtools
    mkdir build
    cd build
    cmake ..
    make
    cd ..
    cp bin/bamtools /usr/local/bin/.

Install novoalign

    - Novoalign is closed-source, commercial software that can be freely used
      by any not-for-profit projects within not-for-profit organizations.
      Visit the 'Downloads' page at http://www.novocraft.com/ 

Install yaha
::
    wget http://faculty.virginia.edu/irahall/support/yaha/YAHA.0.1.72.tar.gz
    tar zxvf YAHA.0.1.72.tar.gz
    cp yaha /usr/local/bin/.

Install bwa
::
    wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.4.tar.bz2
    bunzip2 bwa-0.7.4.tar.bz2
    tar xvf bwa-0.7.4.tar
    cd bwa-0.7.4
    make
    cp bwa /usr/local/bin/.

Clone the Lumpy repository
::

   git clone git://github.com/arq5x/lumpy-sv.git

Navigate into the lumpy directory
::

  cd lumpy-sv

Edit the `defs.local` file in accordance with your configuration.
    - Edit the `GSL_INCLUDE` environment variable.
        * This path depends on where gsl was installed.  The path should
          contain the "gsl_statistics_int.h" file.  Possibilities are:
        * if OSX,   set GSL_INCLUDE=-I/usr/local/include/ -I/usr/local/include/gsl
        * if Linux, set GSL_INCLUDE=-I/usr/include/gsl/
        * if Windows, sorry this is unsupported.
    - Edit the `GSL_LINK` environment variable.
        * This path depends on where gsl was installed.  The path should
          contain the "libgsl.so.0" file.  Possibilities are:
        * if OSX,   set GSL_LINK=-L/usr/local/lib/
        * if Linux, set GSL_LINK=-L/usr/lib/
        * if Windows, sorry this is unsupported.
        
Finally make sure you have the Statistics::Descriptive CPAN package installed for Perl
::
	sudo cpan Statistics::Descriptive

At this point, you should be ready to compile lumpy
::

        make


Now, you can test the paired-end, split-read, and combination paired-end and
split-read  versions of the tools by running the moving to the `test` directory and running the `test.sh` script. Also, this shell script demonstrates how 
to run each of the lumpy version
::

        ./test.sh

If all works well, and you have bedtools installed in your path you should see
the following
::
    Testing lumpy paired-end
    Simulated:1000  Predicted:40    True:40 False:0
    Testing lumpy split-read
    Simulated:1000  Predicted:44    True:44 False:0
    Testing lumpy paired-end and split-read
    Simulated:1000  Predicted:95    True:93 False:0

Usage
=====

General options
::

    -e  Show evidnece for each call

The default output reports the predicted breakpoint.  This option includes the
evidence supporting each call.
::

    -mw minimum weight for a call

Each piece of evidence has a weight, and each possible call has an evidence
set.  The sum of weights in the evidence set must be above this value.
::

    -tt trim threshold

Each prediced breakpoint interval has a probability array associated with it.
The intervals can be trimmed of values that are below some trimming percentile.

Split-read options
::

    -sr 
        bam_file:<file name>,

Position sorted bam file containing the output of a single read split-read
aligner (e.g., YAHA, bwasw) for this sample.
::

        back_distance:<distance>

The distance around the +/- of the split to include in the breakpoint interval.
A distance of 20 will created a breakpoint interval of size 40 centered at the
split.
::

        min_mapping_threshold:<mapping quality>

Minimum mapping quality (reported from the aligner) that a read must have 
to be considered.  A quality of 1 will filter all reads with two or more 
equally good mappings.
::

        weight:<sample weight>

Weight of each piece of evidence from this sample.
::

        id:<sample id>

Sample id.

Paired-end options
::

    -pe 
        bam_file:<file name>,

Position sorted bam file containing the output of a paired-end read aligner
aligner (e.g., bwa) for this sample.
::

        histo_file:<file name>,

Histogram of observed library sizes for the sample.  A script to 
generate this file is located in scripts/run_histo.sh
::

        mean:<value>,

Sample mean library size (can be found using scripts/run_histo.sh)
::

        stdev:<value>,

Sample mean library standard deviation (can be found using scripts/run_histo.sh)
::

        read_length:<length>,

Length of sequenced reads
::

        min_non_overlap:<length>,

Number of base pair positions that must be unique to each end of a read pair.
Some library preps are created with large reads and small library sizes such
that read overlap, in all over cases overlapping reads tends to be a sign of an
error.  We typically set this to read length (pairs cannot overlap).
::

        discordant_z:<z value>,

Number of standard deviations away from the mean to be considered as a normal
library size.
::

        back_distance:<distance>

Distance into the read to add to the breakpoint interval. 
::

        min_mapping_threshold:<mapping quality>

Minimum mapping quality (reported from the aligner) that a read must have 
to be considered.  A quality of 1 will filter all reads with two or more 
equally good mappings.
::

        weight:<sample weight>

Weight of each piece of evidence from this sample.
::

        id:<sample id>

Sample id.



BEDPE (general interface) options
::

    -bedpe 
        bedpe_file:<bedpe file>,

Position sorted bedpe file containing the breakpoint intervals for this sample.
::

        distro_file:<distro_file>,

File containing the values for the breakpoint probability array.
::

        back_distance:<distance>

Distance into the read to add to the breakpoint interval.  
::

        weight:<sample weight>

Weight of each piece of evidence from this sample.
::

        id:<sample id>

Sample id.


Output
======

Tab separated::

	1. chromosome 1
	2. interval 1 start
	3. interval 1 end
	4. chromosome 2
	5. interval 2 start
	6. interval 2 end
	7. id
	8. evidence set size
	9. strand 1
	10. strand 2
	11. type (DELETION = 1, DUPLICATION = 2, INVERSION = 3)
	12. id of samples containing evidence for this breakpoint

Example::

	chr10	2225782	2226073	chr10	2235576	2235865	0x10f504f80	4	+	-	1	ids:1

Test data sets
==============
The `test/test.sh` script executes lumpy against several simulated data sets
and compares the results to the known correct result.  The sample data sets are
not part of the lumpy code base, and can be found at
`http://www.cs.virginia.edu/~rl6sf/lumpy/data.tar.gz`.  This tar ball should be
extracted into the top-level lumpy directory.  The script `test/test.sh` checks
for the the existence of this directory before running lumpy.

Example Single Sample PE and SR Work flow
========================================

Assume that the input files are "sample.1.fq" and "sample.2.fq", and the read length is 150.

Paired end alignment
-----

We prefer novoalign for paired-end reads:
::
    novoalign \
        -d hg19.ndx \
        -o SAM \
        -r Random \
        -i PE 500,50 -e 1 -c 20 \
        -f sample.1.fq sample.2.fq \
        | samtools view -Sb - > sample.pe.bam

However, bwa is another option:
::
    bwa aln hg19.fa sample.1.fq > sample.1.sai
    bwa aln hg19.fa sample.2.fq > sample.2.sai
    bwa sampe hg19.fa \
        sample.1.sai sample.2.sai \
        sample.1.fq sample.2.fq \
        | samtools view -S -b - \
        > sample.pe.bam

Use bamtools to sort since samtools does not correctly flag the resulting bam as coordinate sorted:
::
    bamtools sort -in sample.pe.bam -out sample.pe.sort.bam

Split read alignment
-----

Extract the reads in sample.pe.sort.bam that are either unmapped or have a soft clipped portion of at least 20 base pairs
::
    samtools view sample.pe.sort.bam \
        | scripts/split_unmapped_to_fasta.pl -b 20 \
	> sample.um.fq

Use a split-read aligner on the unmapped/soft clipped reads; we prefer yaha:
::
    # index first
    yaha -g hg19.fa  -L 11
    
    # using 20 threads
    yaha \
        -t 20 \
	-x hg19.X11_01_65525S
	-q sample.um.fq \
	-osh stdout \
	-M 15 \
	-H 2000 \
	-L 11 \
	| samtools view -Sb - \
	> sample.sr.bam

For split reads, bwasw is another option:
::    
    bwa bwasw -H -t 20 hg19.fa sample.um.fq \
        | samtools view -Sb - \
        > sample.sr.bam

Sort the split-read alignments (again, using bamtools):
::
    bamtools sort -in sample.sr.bam -out sample.sr.sort.bam

Run lumpy-sv using paired end reads
-----

Using the paired end mapped reads,  empirically define the paired-end distribution from 10000 proper alignments:
::    
    samtools view sample.pe.sort.bam \
        | scripts/pairend_distro.pl \
        -rl 150 \
        -X 4 \
        -N 10000 \
        -o sample.pe.histo

The above script (scripts/pairend_distro.pl) will display mean and stdev to screen.

To run lumpy with just the paired-end data, We will assume the mean=500 and stdev=50:
::
    ../bin/lumpy \
        -mw 4 \
	-tt 1e-3 \
	-pe \
	bam_file:sample.pe.sort.bam,histo_file:sample.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:1\
	> sample.pe.bedpe

Run lumpy-sv using split-reads reads
-----

We can run lumpy with just the split-read data too:
::    
    ../bin/lumpy \
        -mw 4 \
	-tt 1e-3 \
	-sr \
	bam_file:sample.sr.sort.bam,back_distance:20,weight:1,id:1,min_mapping_threshold:1 \
	> sample.sr.bedpe

Run lumpy-sv using both paired and split reads
-----

Or, we run lumpy with both the paired-end and split-read data:
::
	../bin/lumpy \
		-mw 4 \
		-tt 1e-3 \
		-pe \
		bam_file:sample.pe.sort.bam,histo_file:sample.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:1\
		-sr \
		bam_file:sample.sr.sort.bam,back_distance:20,weight:1,id:1,min_mapping_threshold:1 \
		> sample.pesr.bedpe

Run lumpy-sv using matched samples
-----

Lastly, we can run lumpy with paired-end data from a matched tumor/normal samples
::
	../bin/lumpy \
	        -mw 4 \
	        -tt 1e-3 \
	        -pe \
	        bam_file:tumor.pe.sort.bam,histo_file:tumor.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:1\
	        -pe \
	        bam_file:normal.pe.sort.bam,histo_file:normal.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:1\
	        > tumor_v_normal.pe.bedpe
