#!/usr/bin/env python
#  (c) 2012 - Ryan M. Layer
#  Hall Laboratory
#  Quinlan Laboratory
#  Department of Computer Science
#  Department of Biochemistry and Molecular Genetics
#  Department of Public Health Sciences and Center for Public Health Genomics,
#  University of Virginia
#  rl6sf@virginia.edu

import sys
import numpy as np
from operator import itemgetter
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-r",
    "--read_length",
    type="int",
    dest="read_length",
    help="Read length")

parser.add_option("-X",
    dest="X",
    type="int",
    help="Number of stdevs from mean to extend")

parser.add_option("-N",
    dest="N",
    type="int",
    help="Number to sample")

parser.add_option("-o",
    dest="output_file",
    help="Output file")


(options, args) = parser.parse_args()

if not options.read_length:
    parser.error('Read length not given')

if not options.X:
    parser.error('X not given')

if not options.N:
    parser.error('N not given')

if not options.output_file:
    parser.error('Output file not given')


required = 67
restricted = 384

L = []
c = 0

for l in sys.stdin:
    if c >= options.N:
        break

    A = l.rstrip().split('\t')

    if ( ((int(A[1]) & required) == required) and \
          (int(A[1]) & restricted == 0 ) and \
          (int(A[8]) >= 0) ):
        L.append(float(A[8]))
        c += 1

mean = np.mean(L)
stdev = np.std(L)

start = options.read_length
end = int(mean + options.X*stdev)

H = [0] * (end - start + 1)
s = 0

for i in range(c):
    if (L[i] >= start) and (L[i] <= end):
        j = int(L[i] - start)
        H[j] = H[ int(L[i] - start) ] + 1
        s += 1

f = open(options.output_file, 'w')

for i in range(end - start):
    o = str(i) + "\t" + str(float(H[i])/float(s)) + "\n"
    f.write(o)


f.close()

print 'mean:' + str(mean) + '\tstdev:' + str(stdev)
