#!/usr/bin/python3
import json
import os
import subprocess
import multiprocessing
import sys

# Parse json
with open("/batchx/input/input.json", "r") as inputFile:
    inputJson = inputFile.read()
parsedJson = json.loads(inputJson)

# BX_MEMORY & BX_VCPUS
bxMemory = os.environ['BX_MEMORY']
bxVcpus = os.environ['BX_VCPUS']

# temp directory on BatchX
tmpDir = "/tmp/"

# Define output dir
outputDir = "/batchx/output/"

# outputFolderName
if "outputFolderName" in parsedJson:
    outputFolderName = parsedJson["outputFolderName"]
    outputFolder = outputDir + outputFolderName + "/"
else:
    outputFolder = outputDir + "lumpy-express-results/"

os.mkdir(outputFolder)

# get BAM file
bamFile = parsedJson["bamFile"] 
bamFileString = ",".join(bamFile)

# get bamExtractFiles
if "bamExtractFiles" in parsedJson:
    bamExtractFiles = parsedJson["bamExtractFiles"]
    bamSplitReadsFile = bamExtractFiles["bamSplitReadsFile"]
    bamSplitReadsFileString = ",".join(bamSplitReadsFile)
    bamDiscReadsFile = bamExtractFiles["bamDiscReadsFile"]
    bamDiscReadsFileString = ",".join(bamDiscReadsFile)
    bamExtractFilesString = "-S " + bamSplitReadsFileString + " -D " + bamDiscReadsFileString
else:
    bamExtractFilesString = ""

# get blackListBed
if "blackListBed" in parsedJson:
    blackListBed = parsedJson["blackListBed"]
    blackListBedString = "-x " + blackListBed
else:
    blackListBedString = ""

# get outputProb
if "outputProb" in parsedJson:
    outputProb = parsedJson["outputProb"]
    if outputProb is True:
        outputProbString = "-P"
    elif outputProb is False:
        outputProbString = ""
else:
    outputProbString = ""

# get minSampleWt
if "minSampleWt" in parsedJson:
    minSampleWt = parsedJson["minSampleWt"]
else:
    minSampleWt = 4

# get trimThreshold
if "trimThreshold" in parsedJson:
    trimThreshold = parsedJson["trimThreshold"]
else:
    trimThreshold = 0

# get tempDir
if "tempDir" in parsedJson:
    tempDir = parsedJson["tempDir"]
    if tempDir is True:
        tempDirString = "-T " + tmpDir
    elif tempDir is False:
        tempDirString = ""
else:
    tempDirString = ""

# get keepTempFiles
if "keepTempFiles" in parsedJson:
    keepTempFiles = parsedJson["keepTempFiles"]
    if keepTempFiles is True:
        keepTempFilesString = "-k"
    elif keepTempFiles is False:
        keepTempFilesString = ""
else:
    keepTempFilesString = ""

# get outputPrefix
if "outputPrefix" in parsedJson:
    outputPrefix = parsedJson["outputPrefix"] + ".vcf"
else:
    outputPrefix = "fullBam.bam.vcf"

try:
    # subprocess.check_call ("python --version", shell=True)
    # subprocess.check_call ("lumpy --version", shell=True)
    # lumpyCmd = "lumpyexpress -B " +  bamFileString + " " + bamExtractFilesString + " -o " + outputFolder + outputPrefix + " " + blackListBedString + " " + outputProbString + " -m " + str(minSampleWt) + " -r " + str(trimThreshold) + " " + tempDirString + " " + keepTempFilesString + " -v"
    lumpyCmd = "lumpyexpress -B " +  bamFileString + " -o " + outputFolder + outputPrefix + " -v"
    subprocess.check_call (lumpyCmd, shell=True)
    # subprocess.check_call (f"lumpyexpress -B {bamFileString} {bamExtractFilesString} -o {outputFolder}/{outputPrefix} \
    # {blackListBedString} {outputProbString} -m {minSampleWt} -r {trimThreshold} {tempDirString} {keepTempFilesString} \
    # {lumpyConfigFileString} -v", shell=True)
except subprocess.CalledProcessError as e:
        print(e)
        exit(e.returncode)
except:
    print("Unexpected error:", sys.exc_info()[0])
    raise

# Write output json file
outputJson = {
    'outputFolder': outputFolder
}
with open('/batchx/output/output.json', 'w+') as json_file:
    json.dump(outputJson, json_file)