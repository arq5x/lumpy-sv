#!/usr/local/bin/python
import sys
import subprocess
import os
import numpy as np

if len(sys.argv) < 2:
    print 'usage:' + sys.argv[0] + ' <in bam 1> <in bam 2> <..>'
    exit(1)

for i in range(1,len(sys.argv)):
    bam_file = sys.argv[i]

    p = subprocess.Popen(\
            ['samtools', 'view', '-H', bam_file], \
            stdout=subprocess.PIPE)

    f = open(bam_file + '.genome', 'w')

    for l in p.communicate()[0].split('\n'):
        if '@SQ' in l:
            A = l.rstrip().split('\t')
            name = ''
            size = ''
            for a in A:
                if 'SN' in a:
                    name = a.split(':')[1]
                if 'LN' in a:
                    size = a.split(':')[1]
            f.write(name + '\t' + size + '\n')
    f.close()

    p = subprocess.Popen(\
            'bedtools genomecov ' + \
                ' -ibam ' +  bam_file + \
                ' -g ' + bam_file + '.genome' + \
                ' -bga ' + \
                ' > ' + bam_file + '.coverage',shell=True)
    os.waitpid(p.pid, 0)

for i in range(1,len(sys.argv)):
    bam_file = sys.argv[i]
    genome_file = bam_file + '.genome'
    f = open(genome_file,'r')
    total_len = 0.0
    for l in f:
        total_len += float(l.rstrip().split('\t')[1])
    f.close()

    coverage_file = bam_file + '.coverage'
    f = open(coverage_file, 'r')
    C = []
    W = []
    for l in f:
        a = l.rstrip().split('\t')
        if float(a[3]) > 0:
            C.append(float(a[3])) 
            W.append((float(a[2])-float(a[1]))/total_len)
    min_c = min(C)
    max_c = max(C)
    mean_c = np.average(C,weights=W)
    stdev_c = np.std(C)
    print coverage_file + \
            '\tmin:' + str(min_c) + \
            '\tmax:' + str(max_c) + \
            '\tmean(non-zero):' + str(mean_c) 
    f.close()
