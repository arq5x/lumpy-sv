#!/usr/local/bin/python
import sys
import subprocess
import os
import numpy as np

if len(sys.argv) < 2:
    print 'usage:' + sys.argv[0] + ' <bam 1> <bam 2> <..>'
    exit(1)

order = []


for i in range(1,len(sys.argv)):
    bam_file = sys.argv[i]
    print bam_file

    p = subprocess.Popen(\
            ['samtools', 'view', '-H', bam_file], \
            stdout=subprocess.PIPE)

    for l in p.communicate()[0].split('\n'):
        if '@SQ' in l:
            A = l.rstrip().split('\t')
            name = ''
            for a in A:
                if 'SN' in a:
                    name = a.split(':')[1]
            order.append(name)


    p = subprocess.Popen(\
            'samtools view ' +  bam_file + \
            ' | cut -f3,4 | head -n 100000', \
            shell=True, \
            stdout=subprocess.PIPE)

    broke = False
    curr_chrom_index = -1
    curr_pos = -1
    for l in p.communicate()[0].split('\n'):
        if len(l) > 0 :
            a = l.split('\t')
            chrom = a[0]
            pos = int(a[1])

            if order.index(chrom) > curr_chrom_index:
                curr_chrom_index = order.index(chrom)
                curr_pos = -1
            elif order.index(chrom) < curr_chrom_index:
                print 'out of order:\t' + l + '\toccurred after\t' + \
                        order[curr_chrom_index] + '\t' + str(curr_pos)
                broke = True
                break

            if pos > curr_pos:
                curr_pos = pos
            elif pos < curr_pos:
                print 'out of order:\t' + l + '\toccurred after\t' + \
                        order[curr_chrom_index] + '\t' + str(curr_pos)
                broke = True
                break
    if not broke:
        print "in order"


