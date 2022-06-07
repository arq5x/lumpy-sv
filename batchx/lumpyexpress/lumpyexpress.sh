#!/usr/bin/env bash
set -eo pipefail

# Input json
input=/batchx/input/input.json

# Input BAMs
tmpBams=/tmp/tmpBams/
mkdir $tmpBams

# BAM Files
bamFilesString=''
for bamFile in $(cat $input | jq -r .bamBaiFiles.bamFile[]);
do
    ln -s $bamFile "$tmpBams"
    bamFile="$tmpBams$(basename -- "$bamFile")"
	if [ "${#bamFilesString}" -gt 0 ]; then
		bamFilesString="$bamFilesString,"
	fi
    bamFilesString="$bamFilesString$bamFile"
done
numberOfBamFiles=$(cat $input | jq -r .bamBaiFiles.bamFile[] | wc -l)

# BAI Files
baiFiles=$(cat $input | jq -r .bamBaiFiles.baiFile[]?)
if [ "$baiFiles" != "" ]; then
    for baiFile in "$baiFiles";
    do
        ln -s $baiFile "$tmpBams"
        baiFile="$tmpBams$(basename -- "$baiFile")"
    done
fi

# get BAM Extract files
bamExtractFiles=$(cat $input | jq -r .bamExtractFiles)
if [ "$bamExtractFiles" != 'null' ]; then
    # bamSplitReadsFile
    bamSplitReadsFileString=''
    for bamSplitFile in $(cat $input | jq -r .bamExtractFiles.bamSplitReadsFile[]);
    do
	    if [ "${#bamSplitReadsFileString}" -gt 0 ]; then
		    bamSplitReadsFileString="$bamSplitReadsFileString,"
	    fi
        bamSplitReadsFileString="$bamSplitReadsFileString$bamSplitFile"
    done
    numberOfBamSplitReadsFiles=$(cat $input | jq -r .bamExtractFiles.bamSplitReadsFile[] | wc -l)
    # bamDiscReadsFile
    bamDiscReadsFileString=''
    for bamDiscFile in $(cat $input | jq -r .bamExtractFiles.bamDiscReadsFile[]);
    do
	    if [ "${#bamDiscReadsFileString}" -gt 0 ]; then
		    bamDiscReadsFileString="$bamDiscReadsFileString,"
	    fi
        bamDiscReadsFileString="$bamDiscReadsFileString$bamDiscFile"
    done
    numberOfBamDiscReadsFiles=$(cat $input | jq -r .bamExtractFiles.bamDiscReadsFile[] | wc -l)
    if [ $numberOfBamSplitReadsFiles != $numberOfBamDiscReadsFiles ]; then
        echo "The number of 'bamSplitReadsFile' and 'bamDiscReadsFile' files provided is not the same. Each splitter reads BAM file should have a corresponding discordant reads BAM file and vice-versa. Please re-check the 'bamSplitReadsFile' and 'bamDiscReadsFile' files, and provide them again."
        exit 1
    fi
    if [ $numberOfBamFiles != $numberOfBamSplitReadsFiles ]; then
        echo "The number of 'bamSplitReadsFile' and 'bamFile' files provided is not the same. Each BAM file should have a splitter reads and a discordant reads BAM file. Please re-check the 'bamSplitReadsFile' and 'bamDiscReadsFile' files, and provide them again."
        exit 1
    fi
    bamExtractFilesCmd="-S $bamSplitReadsFileString -D $bamDiscReadsFileString"
fi

# get blackListBed
blackListBed=$(cat $input | jq -r .blackListBed)
if [ $blackListBed != 'null' ]; then
    blackListBedCmd="-x $blackListBed" 
fi

# get outputProb
outputProb=$(cat $input | jq -r .outputProb)
if [ $outputProb != 'null' ] && [ $outputProb = true ]; then
    outputProbCmd='-P'
fi

# get minSampleWt
minSampleWt=$(cat $input | jq -r .minSampleWt)
if [ $minSampleWt = 'null' ]; then
    minSampleWt=4
fi

# get trimThreshold
trimThreshold=$(cat $input | jq -r .trimThreshold)
if [ $trimThreshold = 'null' ]; then
    trimThreshold=0
fi

# get outputPrefix

outputPrefix=$(cat $input | jq -r .outputPrefix)
if [ $outputPrefix = 'null' ]; then
    outputPrefix="fullBam.bam"
fi

# Output dir
outputFolder=/batchx/output/lumpy-express-results/
mkdir $outputFolder

# vcfFile
vcfFile="$outputFolder$outputPrefix.vcf"

# get keepTempFiles
keepTempFiles=$(cat $input | jq -r .keepTempFiles)
if [ $keepTempFiles != 'null' ] && [ $keepTempFiles = true ]; then
    # temp directory on BatchX where output temporary files would go if keepTempFiles is True
    tempFiles=/batchx/output/lumpy-express-results/tempFiles
    mkdir $tempFiles
    outputTmpDirCmd="-T $tempFiles"
    keepTempFilesCmd='-k'
fi

# Run lumpyexpress
echo "Running lumpyexpress..."
lumpyexpress -B $bamFilesString $bamExtractFilesCmd -o $vcfFile $blackListBedCmd $outputProbCmd \
 -m $minSampleWt -r $trimThreshold $keepTempFilesCmd $outputTmpDirCmd -v
echo "Finished running lumpyexpress."

# Output json
if [ $keepTempFiles != 'null' ] && [ $keepTempFiles = true ]; then
    echo "Outputting temporary files..."
    echo "{\"vcfFile\":\"$vcfFile\",\"tempFiles\":\"$tempFiles\"}" > /batchx/output/output.json    
else
    echo "{\"vcfFile\":\"$vcfFile\"}" > /batchx/output/output.json
fi

