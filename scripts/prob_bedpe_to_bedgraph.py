#!/usr/bin/env python
import sys
import numpy as np

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-i",
    "--bedpe",
    dest="bedpe",
    help="BEDPE file")

parser.add_option("-n",
    "--name",
    default="LUMPY BedGraph",
    dest="name",
    help="Name")


(options, args) = parser.parse_args()

if not options.bedpe:
    parser.error('BEDPE file not given')

f = open(options.bedpe,'r')

print 'track type=bedGraph name="' + options.name + '"' 

for line in f:
    if line[0] == '#':
        continue

    v = line.rstrip().split('\t')
    info = dict()
    for s in v[12].split(';'):
        s_v = s.split('=')
        if len(s_v) == 2:
            info[s_v[0]] = s_v[1]

    L=[float(x) for x in info['PRPOS'].split(',')] 
    R=[float(x) for x in info['PREND'].split(',')] 

    l_chr = v[0]
    l_start = int(v[1])

    r_chr = v[3]
    r_start = int(v[4])

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
