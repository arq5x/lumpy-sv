#!/usr/local/bin/python
import sys
import subprocess
import os
import numpy as np

if len(sys.argv) < 3:
    print 'usage:' + sys.argv[0] + ' <max> <out file> <in bam 1> <in bam 2> <..>'
    exit(1)

max_c = int(sys.argv[1])
out_file = sys.argv[2]

o = open('.exclude.tmp','w')

for i in range(3,len(sys.argv)):
    bam_file = sys.argv[i]
    coverage_file = bam_file + '.coverage'

    f = open(coverage_file, 'r')
    for l in f:
        a = l.rstrip().split('\t')
        if (int(a[3]) > max_c):
            o.write(l)
    f.close()
o.close()
            
p = subprocess.Popen(\
        'cat .exclude.tmp | ' \
        'sort -S 20G -k1,1 -k2,2n | ' \
        'bedtools merge -i stdin -d 50 ' \
        '> ' + out_file, shell=True)
os.waitpid(p.pid, 0)
os.unlink('.exclude.tmp')
