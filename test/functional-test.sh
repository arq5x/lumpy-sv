#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -e

make

if [[ ! -e data/pe.pos_sorted.bam ]]; then

	wget http://layerlab.org/lumpy/data.tar.gz
	tar xzvf data.tar.gz
fi

_fn_bp_pe_sr() {
	set -eo nounset

MIN_MAP_T=$1
BACK=$2
TT=$3
READ_LENGTH=$4
PE_BAM=$5
SR_BAM=$6
E=$7
Z=4
WEIGHT=4

OUT_DIR=`dirname $PE_BAM`
OUT_FILE=`basename $PE_BAM .bam`

MEAN_STDEV=`samtools view $PE_BAM \
				| ./scripts/pairend_distro.py \
					-r $READ_LENGTH \
					-X $Z \
					-N 10000 \
					-o $OUT_DIR/$OUT_FILE.histo`


MEAN=`echo $MEAN_STDEV | cut -d " " -f1 | cut -d ":" -f2`
STDEV=`echo $MEAN_STDEV | cut -d " " -f2 | cut -d ":" -f2`

./bin/lumpy \
    -b \
	-mw $WEIGHT \
    -tt $TT \
    -sr \
    bam_file:$SR_BAM,back_distance:$2,weight:1,id:1,min_mapping_threshold:$MIN_MAP_T \
    -pe \
    bam_file:$PE_BAM,histo_file:$OUT_DIR/$OUT_FILE.histo,mean:$MEAN,stdev:$STDEV,read_length:$READ_LENGTH,min_non_overlap:$READ_LENGTH,discordant_z:4,back_distance:$BACK,weight:1,id:1,min_mapping_threshold:$MIN_MAP_T\
    > $OUT_DIR/$OUT_FILE.pesr.bedpe

if [ `command -v bedtools` ]
then
	EX=`cat $E | wc -l`

	O=`cat $OUT_DIR/$OUT_FILE.pesr.bedpe | wc -l`

	C=`bedtools pairtopair \
		-a $E \
		-b $OUT_DIR/$OUT_FILE.pesr.bedpe \
		-type both -is -slop 3 \
		| cut -f7 | sort -u | wc -l`

	I=`bedtools pairtopair \
		-b $E \
		-a $OUT_DIR/$OUT_FILE.pesr.bedpe \
		-type notboth -is -slop 3 \
		| cut -f 7 | sort -u | wc -l`

	echo -e "Simulated:$EX\tPredicted:$O\tTrue:$C\tFalse:$I"
else
	echo "Install bedtools to compare result to known deletions"
fi

}

fn_bp_pe_sr() {
	_fn_bp_pe_sr 1 20 1e-3 150 data/pe.pos_sorted.bam \
			    data/sr.pos_sorted.bam \
				data/sim.bedpe
}

run lumpy_test fn_bp_pe_sr
assert_exit_code 0
cat $STDOUT_FILE
