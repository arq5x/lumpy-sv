#!/bin/bash


if [ ! -d "../data/" ]
then
	echo "../data/ directory does not exist"
	echo "please download from http://www.cs.virginia.edu/~rl6sf/lumpy/data.tar.gz"
	exit
fi


echo "Testing lumpy paired-end"
./bp_pe_pl.sh 1 20 1e-3 150 \
	../data/pe.pos_sorted.bam \
	../data/sim.bedpe

echo "Testing lumpy split-read"
./bp_sr.sh 1 20 1e-3 150 \
	../data/sr.pos_sorted.bam \
	../data/sim.bedpe

echo "Testing lumpy paired-end and split-read"
./bp_pesr_pl.sh 1 20 1e-3 150 \
	../data/pe.pos_sorted.bam \
	../data/sr.pos_sorted.bam \
	../data/sim.bedpe
