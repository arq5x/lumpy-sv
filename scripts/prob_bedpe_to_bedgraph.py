#!/usr/bin/env python
import sys
import numpy as np

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-b",
    "--bedpe_file",
    dest="bedpe_file",
    help="BEDPE file")

parser.add_option("-n",
    "--name",
    default="LUMPY BedGraph",
    dest="name",
    help="Name")


(options, args) = parser.parse_args()

if not options.bedpe_file:
    parser.error('BEDPE file not given')

f = open(options.bedpe_file,'r')

print 'track type=bedGraph name="' + options.name + '"' 

for l in f:
    A = l.rstrip().split('\t')
    L=[float(x) for x in A[15][3:].split(',')] 
    R=[float(x) for x in A[16][3:].split(',')] 

    l_chr = A[0]
    l_start = int(A[1])

    r_chr = A[3]
    r_start = int(A[4])

    c = 0
    for p in L:
        print '\t'.join( [l_chr,
                          str(l_start + c),
                          str(l_start + c + 1),
                          str(p)])
        c+=1

    c = 0
    for p in R:
        print '\t'.join( [r_chr,
                          str(r_start + c),
                          str(r_start + c + 1),
                          str(p)])
        c+=1
    

f.close()
