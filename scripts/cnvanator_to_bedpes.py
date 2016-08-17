#!/usr/bin/env python
import sys
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

parser.add_option("--cnvkit",
        action="store_true",
        default=False,
        help="input is .cns file from cnvkit")

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

cnvkit = options.cnvkit

i = 1
for l in f:
    A = l.rstrip().split('\t')

    if cnvkit:
        # skip header.
        if i == 1 and A[1] == 'start': continue

        # see: https://www.biostars.org/p/178685/#203630
        # require a log2 change > +/- 0.2
        FC = float(A[4])
        if abs(FC) < 0.2: continue

        ev = 'DUPLICATION'

        # check the log2 change.
        if FC < 0:
            ev = 'DELETION'
        call = "%s:%s-%s" % (A[0], A[1], A[2])
        bedpe = interval_to_bedpe(options.breakpoint_size, call, ev, i)
    else:
        ev = A[0].upper()
        bedpe = interval_to_bedpe(options.breakpoint_size, A[1], ev, i)

    assert ev in ("DUPLICATION", "DELETION"), ev

    if ev == "DUPLICATION":
        dup_f.write(bedpe + '\n')
    else:
        del_f.write(bedpe + '\n')
    i += 1

f.close()
del_f.close()
dup_f.close()
