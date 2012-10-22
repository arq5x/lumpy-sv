#!/bin/bash

#  run_histo.sh
#  (c) 2012 - Ryan M. Layer
#  Hall Laboratory
#  Quinlan Laboratory
#  Department of Computer Science
#  Department of Biochemistry and Molecular Genetics
#  Department of Public Health Sciences and Center for Public Health Genomics,
#  University of Virginia
#  rl6sf@virginia.edu
# 
# Licenced under the GNU General Public License 2.0 license.



if [ -z $4 ]
then
	echo "usage:$0 <bam file> <out file> <Z> <read lengh>"
	exit
fi

INF=$1
OUTF=$2
CZ=$3
RL=$4

SCRIPT_DIR=`dirname $0`

OUT=`samtools view $INF\
	| $SCRIPT_DIR/pairend_distro.pl \
	-rl 150 \
	-X $CZ \
	-N 100000 \
	-o $OUTF`

MEAN=`echo $OUT | cut -d" " -f1 | cut -d":" -f2`
STDEV=`echo $OUT | cut -d" " -f2 | cut -d":" -f2`

echo $MEAN
echo $STDEV
