##!/bin/bash

if [ -z $6 ]
then
	echo "$0 <min_mapping_threshold> <back_distance> <tt> <read len> <pe bam> <sr bam> "
	exit
fi

MIN_MAP_T=$1
BACK=$2
TT=$3
READ_LENGTH=$4
PE_BAM=$5
SR_BAM=$6
Z=4
WEIGHT=4

OUT_FILE=`basename $PE_BAM .bam`

MEAN_STDEV=`samtools view $PE_BAM \
| ../../scripts/pairend_distro.py \
    -r $READ_LENGTH \
    -X $Z \
    -N 10000 \
    -o $OUT_FILE.histo`


MEAN=`echo $MEAN_STDEV | cut -d " " -f1 | cut -d ":" -f2`
STDEV=`echo $MEAN_STDEV | cut -d " " -f2 | cut -d ":" -f2`

../../bin/lumpy \
    -b \
	-mw $WEIGHT \
    -tt $TT \
    -sr \
    bam_file:$SR_BAM,back_distance:$2,weight:1,id:1,min_mapping_threshold:$MIN_MAP_T \
    -pe \
    bam_file:$PE_BAM,histo_file:$OUT_FILE.histo,mean:$MEAN,stdev:$STDEV,read_length:$READ_LENGTH,min_non_overlap:$READ_LENGTH,discordant_z:4,back_distance:$BACK,weight:1,id:1,min_mapping_threshold:$MIN_MAP_T
