##!/bin/bash

if [ -z $5 ]
then
	echo "$0 <min_mapping_threshold> <back_distance> <tt> <read len> <rl bam> <e>"
	exit
fi

MIN_MAP_T=$1
BACK=$2
TT=$3
READ_LENGTH=$4
SR_BAM=$5
E=$6
Z=4
WEIGHT=4

OUT_DIR=`dirname $SR_BAM`
OUT_FILE=`basename $SR_BAM .bam`

../bin/lumpy \
	-mw $WEIGHT \
    -tt $TT \
    -sr \
    bam_file:$SR_BAM,back_distance:$2,weight:1,id:1,min_mapping_threshold:$1 \
    > $OUT_DIR/$OUT_FILE.sr.bedpe

if [ `command -v bedtools` ]
then
	EX=`cat $E | wc -l`

	O=`cat $OUT_DIR/$OUT_FILE.sr.bedpe | wc -l`

	C=`bedtools pairtopair \
		-a $E \
		-b $OUT_DIR/$OUT_FILE.sr.bedpe \
		-type both -is -slop 3 \
		| cut -f7 | sort -u | wc -l`

	I=`bedtools pairtopair \
		-b $E \
		-a $OUT_DIR/$OUT_FILE.sr.bedpe \
		-type notboth -is -slop 3 \
		| cut -f 7 | sort -u | wc -l`

	echo -e "Simulated:$EX\tPredicted:$O\tTrue:$C\tFalse:$I"
else
	echo "Install bedtools to compare result to known deletions"
fi

