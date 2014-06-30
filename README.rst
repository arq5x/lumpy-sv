Overview
========

This repository contains the source code for the Lumpy structural variation
detection framework developed in the Hall and Quinlan laboratories at the
University of Virginia.

Installation
============
Lumpy installation does not require any external softwart. We do recommend
installing samtools, bedtools, bamtools, novoalign (or bwa), and yaha.  Below
is a step-by-step tutorial for how to install and use Lumpy. This guide assumes
that `/usr/local/bin` is writable and in your path.  If either is not true, use
another directory that is both writable and in your path, or contact your
administrator.  If you have questions, email me.

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

    wget http://faculty.virginia.edu/irahall/support/yaha/YAHA.0.1.79.tar.gz
    tar zxvf YAHA.0.1.79.tar.gz
    cp yaha /usr/local/bin/.

Install bwa
::

    wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.5a.tar.bz2
    bunzip2 bwa-0.7.5a.tar.bz2
    tar xvf bwa-0.7.5a.tar
    cd bwa-0.7.5a
    make
    cp bwa /usr/local/bin/.

Clone the Lumpy repository
::

   git clone git://github.com/arq5x/lumpy-sv.git

Navigate into the lumpy directory
::

  cd lumpy-sv


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
	549	0
	chr10	1000000
	Simulated:    1000	Predicted:      40	True:      40	False:       0
	Testing lumpy split-read
	chr10	1000000
	Simulated:    1000	Predicted:      44	True:      43	False:       1
	Testing lumpy paired-end and split-read
	549	0
	chr10	1000000
	Simulated:    1000	Predicted:      95	True:      92	False:       1

Usage
=====

General options
::

    -e  

The default output reports the predicted breakpoint.  This option includes the
evidence supporting each call.
::

    -mw minimum weight for a call

Each piece of evidence has a weight, and each possible call has an evidence
set.  The sum of weights in the evidence set must be above this value.
::

    -tt trim threshold

Each predicted breakpoint interval has a probability array associated with it.
The intervals can be trimmed of values that are below some trimming percentile.
NOTE: We recommend "-tt 0.0" (no trimming) since LUMPY now reports both the 95%
confidence interval and the most probable single base for each breakpoint.
::

    -P 

Print the breakpoint probability array.
::

    -x excluded regions bed file

Regions of the genome may be excluded from consideration by included them in
bed file format.  Any alignment that overlaps any of the regions will be
ignored.  This is particularly useful when a sample has regions with either too
very low or very high coverage due to biases in sequencing or alignment.  See
below for help creating this file.
::

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
generate this file is located in scripts/pairend_distro.py
::

        mean:<value>,

Sample mean library size (can be found using scripts/pairend_distro.py)
::

        stdev:<value>,

Sample mean library standard deviation (can be found using scripts/pairend_distro.py)
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
	8. evidence set score
	9. strand 1
	10. strand 2
	11. type 
	12. id of samples containing evidence for this breakpoint
        13. strand configurations observed in the evidence set
        14. point within the two breakpoint with the maximum probability
        15. segmetn of each breakpoint that contains 95% of the probability

Example::

        chr1	547154	547462	chr1	547265	547569	1	0.00254453	+	-	TYPE:DELETION	IDS:10,6	STRANDS:+-,6	MAX:chr1:547175;chr1:547569	95:chr1:547169-547225;chr1:547266-547569

Test data sets
==============
The `test/test.sh` script executes lumpy against several simulated data sets
and compares the results to the known correct result.  The sample data sets are
not part of the lumpy code base, and can be found at
`http://www.cs.virginia.edu/~rl6sf/lumpy/data.tar.gz`.  This tar ball should be
extracted into the top-level lumpy directory.  The script `test/test.sh` checks
for the the existence of this directory before running lumpy.

Example Work flow
========================================

Assume that the input files are "sample.1.fq" and "sample.2.fq", and the read length is 150.

LUMPY is designed to consider both paired-end and split-read alignments, and can also consider each independently.  There are two strategies for extracting constructing a split-read bam file that are fully explained below.  One option is to first align a fastq file with a paired-end aligned (novoalign or bwa), extract candidate split reads from those alignments, then realign those candidate reads using a split-read aligner (yaha or bwasw).  If you are starting with an aligned file (e.g., a bam file), this is probably your best option since it does not require full realignment.  Another option is to align using bwa-mem, which will produce both paired-end alignments and split-read alignments in a single pass.  Then, you can split this file into a paired-end file and a split-read file.  This is probably the best option when starting from a fastq file.

Paired-end alignment
-----

Both novoalign and bwa are options for paired-end alignment:
::

    novoalign \
        -d hg19.ndx \
        -o SAM \
        -r Random \
        -i PE 500,50 -e 1 -c 20 \
        -f sample.1.fq sample.2.fq \
        | samtools view -Sb - > sample.pe.bam

    bwa aln hg19.fa sample.1.fq > sample.1.sai
    bwa aln hg19.fa sample.2.fq > sample.2.sai
    bwa sampe hg19.fa \
        sample.1.sai sample.2.sai \
        sample.1.fq sample.2.fq \
        | samtools view -S -b - \
        > sample.pe.bam

Use bamtools or a recent version of samtools (0.1.19) to sort.  NOTE: the resulting bam file must have the coordinate sort flag set (i.e., @HD VN:1.3  SO:coordinate).
::

    bamtools sort -in sample.pe.bam -out sample.pe.sort.bam

    samtools sort sample.pe.bam sample.pe.sort

Split read alignment
-----

From the paired end aligned bam file sample.pe.sort.bam, you can extract the reads that are either unmapped or have a soft clipped portion of at least 20 base pairs
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

Sort the split-read alignments (again, using bamtools or samtools):
::

    bamtools sort -in sample.sr.bam -out sample.sr.sort.bam

    samtools sort sample.sr.bam sample.sr.sort

Paired-end and split-read alignment using bwa-mem
-----

bwa-mem produces a single bam file with both paired-end alignments and split-read alignments
::

    bwa mem hg19.fa sample.1.fq sample.2.fq -M \
        | samtools view -S -b - \
        > sample.pesr.bam

extract the disordant paired-end alignments.
::

    samtools view -u -F 0x0002 sample.pesr.bam  \
        |  samtools view -u -F 0x0100 - \
        | samtools view -u -F 0x0004 - \
        | samtools view -u -F 0x0008 - \
        | samtools view -b -F 0x0400 - \
        > sample.discordant.pe.bam

extract the split-read alignments
::

    samtools view -h sample.pesr.bam \
        | scripts/extractSplitReads_BwaMem -i stdin \
        | samtools view -Sb - \
        > sample.sr.bam

Sort both alignments (again, using bamtools or samtools):
::

    bamtools sort -in sample.discordant.pe.bam -out sample.discordant.pe.sort.bam
    bamtools sort -in sample.sr.bam -out sample.sr.sort.bam

    samtools sort sample.discordant.pe.bam sample.discordant.pe.sort
    samtools sort sample.sr.bam sample.sr.sort


Run lumpy-sv using paired end reads
-----

Using the paired end mapped reads,  empirically define the paired-end distribution from 10000 proper alignments.  It is common practice to skip the first million reads.
::   

    samtools view sample.pesr.bam \
        | tail -n+100000 \
        | scripts/pairend_distro.py \
        -rl 150 \
        -X 4 \
        -N 10000 \
        -o sample.pe.histo

The above script (scripts/pairend_distro.py) will display mean and stdev to screen.

To run lumpy with just the paired-end data, We will assume the mean=500 and stdev=50:
::

    ../bin/lumpy \
        -mw 4 \
	-tt 0.0 \
	-pe \
	bam_file:sample.discordant.pe.sort.bam,histo_file:sample.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:20\
	> sample.pe.bedpe

Run lumpy-sv using split-reads reads
-----

We can run lumpy with just the split-read data too:
::    

    ../bin/lumpy \
        -mw 4 \
	-tt 0.0 \
	-sr \
	bam_file:sample.sr.sort.bam,back_distance:20,weight:1,id:2,min_mapping_threshold:20 \
	> sample.sr.bedpe

Run lumpy-sv using both paired and split reads
-----

Or, we run lumpy with both the paired-end and split-read data:
::

	../bin/lumpy \
		-mw 4 \
		-tt 0.0 \
		-pe \
		bam_file:sample.discordant.pe.sort.bam,histo_file:sample.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:20\
		-sr \
		bam_file:sample.sr.sort.bam,back_distance:20,weight:1,id:2,min_mapping_threshold:20 \
		> sample.pesr.bedpe

Run lumpy-sv using matched samples
-----

We can run lumpy with paired-end data from a matched tumor/normal samples
::

	../bin/lumpy \
	        -mw 4 \
	        -tt 0.0 \
	        -pe \
	        bam_file:tumor.pe.sort.bam,histo_file:tumor.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:1\
	        -pe \
	        bam_file:normal.pe.sort.bam,histo_file:normal.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:2,min_mapping_threshold:1\
	        > tumor_v_normal.pe.bedpe

Run lumpy-sv with regions of very high coverage excluded
-----
We can direct lumpy to ignore certain regions by using the exclude region
option.  In this example we find and then exclude regions that have very high
coverage.  First we use the get_coverages.py script to find the min, max, and
mean coverages of the the sr and pe bam files, and to create coverage profiles
for both files.
::

        python ../scripts/get_coverages.py \
                sample.pe.sort.bam \
                sample.sr.sort.bam

        sample.pe.sort.bam.coverage  min:1   max:14  mean(non-zero):2.35557521272
        sample.sr.sort.bam.coverage  min:1   max:7   mean(non-zero):1.08945936729

From this output, we will choose to exclude regions that have more than 10x
coverage.  To create the exclude file we will use the get_exclude_regions.py
script to create the exclude.bed file
::

        python ../scripts/get_exclude_regions.py \
                10 \
                exclude.bed \
                sample.pe.sort.bam \
                sample.sr.sort.bam
        
We now rerun lumpy with the exclude (-x) option 
::

	../bin/lumpy \
		-mw 4 \
		-tt 0.0 \
                -x exclude.bed \
		-pe \
		bam_file:sample.pe.sort.bam,histo_file:sample.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:1\
		-sr \
		bam_file:sample.sr.sort.bam,back_distance:20,weight:1,id:2,min_mapping_threshold:1 \
		> sample.pesr.exclude.bedpe

Troubleshooting
============
All of the bam files that lumpy processes must be position sorted.  To check if your bams are sorted correctly, use the check_sorting.py script
::

        python ../scripts/check_sorting.py \
                pe.pos_sorted.bam \
                sr.pos_sorted.bam \
                pe.name_sorted.bam
        pe.pos_sorted.bam
        in order
        sr.pos_sorted.bam
        in order
        pe.name_sorted.bam
        out of order:   chr10   102292476   occurred after   chr10   102292893
