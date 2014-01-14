#!/usr/bin/env python
import sys
import numpy as np
from operator import itemgetter
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-b",
    "--bedpe_file",
    dest="bedpe_file",
    help="BEDPE file")

parser.add_option("-g",
    "--genome_file",
    dest="genome_file",
    help="Genome file")


(options, args) = parser.parse_args()

if not options.bedpe_file:
    parser.error('BEDPE file not given')

f = open(options.bedpe_file,'r')

B = {}

for l in f:
    A = l.rstrip().split('\t')

    if not A[0] in B:
        B[A[0]] = []

    B[A[0]].append([int(A[1]), l.rstrip()])

f.close()

order=[]

if options.genome_file:
    f = open(options.genome_file,'r')

    for l in f:
        A = l.rstrip().split('\t')
        order.append(A[0])

    f.close()
else:
    order = sorted(B.keys())



for c in order:
    if c in B:
        for l in sorted(B[c], key=itemgetter(0)):
            print l[1]
