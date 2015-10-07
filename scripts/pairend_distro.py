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

# some constants for sam/bam field ids
SAM_FLAG = 1
SAM_REFNAME = 2
SAM_MATE_REFNAME = 6
SAM_ISIZE = 8

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

parser.add_option("-m",
    dest="mads",
    type="int",
    default=10,
    help="Outlier cutoff in # of median absolute deviations (unscaled, upper only)")


def unscaled_upper_mad(xs):
    """Return a tuple consisting of the median of xs followed by the
    unscaled median absolute deviation of the values in xs that lie
    above the median.
    """
    med = np.median(xs)
    if len(xs[xs > med]):
        umad = np.median(xs[xs > med] - med)
    else:
        # When lengths are all the same, return a safe default.
        # This is common with simulated data used for testing.
        umad = 1
    return med, umad


(options, args) = parser.parse_args()

if not options.read_length:
    parser.error('Read length not given')

if not options.X:
    parser.error('X not given')

if not options.N:
    parser.error('N not given')

if not options.output_file:
    parser.error('Output file not given')


required = 97
restricted = 3484
flag_mask = required | restricted

L = []
c = 0

for l in sys.stdin:
    if c >= options.N:
        break

    A = l.rstrip().split('\t')
    flag = int(A[SAM_FLAG])
    refname = A[SAM_REFNAME]
    mate_refname = A[SAM_MATE_REFNAME]
    isize = int(A[SAM_ISIZE])

    want = mate_refname == "=" and flag & flag_mask == required and isize >= 0
    if want:
        c += 1
        L.append(isize)

# Remove outliers
L = np.array(L)
L.sort()
med, umad = unscaled_upper_mad(L)
upper_cutoff = med + options.mads * umad
L = L[L < upper_cutoff]
new_len = len(L)
removed = c - new_len
sys.stderr.write("Removed %d outliers with isize >= %d\n" %
    (removed, upper_cutoff))
c = new_len

mean = np.mean(L)
stdev = np.std(L)

start = options.read_length
end = int(mean + options.X*stdev)

H = [0] * (end - start + 1)
s = 0

for x in L:
    if (x >= start) and (x <= end):
        j = int(x - start)
        H[j] = H[ int(x - start) ] + 1
        s += 1

f = open(options.output_file, 'w')

for i in range(end - start):
    o = str(i) + "\t" + str(float(H[i])/float(s)) + "\n"
    f.write(o)


f.close()

print('mean:' + str(mean) + '\tstdev:' + str(stdev))
