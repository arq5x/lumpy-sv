../../bin/lumpy -P \
-msw 4 \
-tt 0 \
-pe bam_file:test_sub.filtered_disc.bam,histo_file:test.sample1.lib1.x4.histo,mean:438.420989806,stdev:92.6149432404,read_length:100,min_non_overlap:100,discordant_z:5,back_distance:10,weight:1,id:test,min_mapping_threshold:20,read_group:test \
-sr bam_file:test.filtered_split.bam,back_distance:10,min_mapping_threshold:20,weight:1,id:test,min_clip:20 
