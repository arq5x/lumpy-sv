#!/bin/bash

../../bin/lumpy \
    -mw 4 \
    -tt 0.0 \
    -pe \
    bam_file:AL87.discordant.sort.bam,histo_file:AL87.histo,mean:429,stdev:84,read_length:83,min_non_overlap:83,discordant_z:4,back_distance:1,weight:1,id:1,min_mapping_threshold:20 \
    -sr \
    bam_file:AL87.sr.sort.bam,back_distance:1,weight:1,id:2,min_mapping_threshold:20 
