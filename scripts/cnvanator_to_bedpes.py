#!/usr/bin/env python
import sys
import numpy as np
from optparse import OptionParser

def interval_to_bedpe(size, call, sv_type, i):
    half = size/2
    chrom,interval = call.split(':')
    start,end = [int(x) for x in interval.split('-')]
    length = end-start

    bedpe = '\t'.join([chrom,
                      str( max(1,start-half)),
                      str( start+half),
                      chrom,
                      str( max(1,end-half)),
                      str( end+half),
                      str(i),
                      str(length),
                      '+',
                      '+',
                      'TYPE:' + sv_type])

    return bedpe

parser = OptionParser()

parser.add_option("-c",
    "--cnv_calls",
    dest="cnv_calls",
    help="Output file from CNVanator")

parser.add_option(
    "--del_o",
    dest="del_o",
    help="Deletion output bedpe file name")

parser.add_option(
    "--dup_o",
    dest="dup_o",
    help="Duplication output bedpe file name")

parser.add_option("-b",
    "--breakpoint_size",
    dest="breakpoint_size",
    type="int",
    help="The total size of the resulting breakpoint, " + \
            "centered on the call edge")

(options, args) = parser.parse_args()

if not options.cnv_calls:
    parser.error('CNVanator calls not given')

if not options.del_o:
    parser.error('Deletion output file not given')

if not options.dup_o:
    parser.error('Duplication output file not given')

if not options.breakpoint_size:
    parser.error('Breakpoint size not given')


f = open(options.cnv_calls,'r')

del_f = open(options.del_o,'w')
dup_f = open(options.dup_o,'w')

i = 1
for l in f:
    A = l.rstrip().split('\t')

    if A[0] == 'deletion':
        bedpe = interval_to_bedpe(options.breakpoint_size, A[1],'DELETION',i)
        del_f.write(bedpe + '\n')
    elif A[0] == 'duplication':
        bedpe = interval_to_bedpe(options.breakpoint_size, A[1],'DUPLICATION',i)
        dup_f.write(bedpe + '\n')

    i += 1


f.close()
del_f.close()
dup_f.close()
