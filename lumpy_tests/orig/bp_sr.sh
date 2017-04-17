##!/bin/bash

if [ -z $4 ]
then
	echo "$0 <min_mapping_threshold> <back_distance> <tt> <read len> <rl bam>"
	exit
fi

MIN_MAP_T=$1
BACK=$2
TT=$3
READ_LENGTH=$4
SR_BAM=$5
Z=4
WEIGHT=4

OUT_FILE=`basename $SR_BAM .bam`

../../bin/lumpy \
    -b \
    -mw $WEIGHT \
    -tt $TT \
    -sr \
    bam_file:$SR_BAM,back_distance:$2,weight:1,id:1,min_mapping_threshold:$1 
