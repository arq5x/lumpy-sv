#!/usr/bin/env python
import sys
import numpy as np
import glob
from operator import add
from optparse import OptionParser

class breakpoint:
    chr_l = ''
    start_l = 0
    end_l = 0
    p_l = []

    chr_l = ''
    start_l = 0
    end_l = 0
    p_l = []

    sv_type = ''

    evidence_size = 0


    def __init__(self, line):
        A = l.rstrip().split('\t')
        self.chr_l = A[0]
        self.start_l = int(A[1])
        self.end_l = int(A[2])

        self.chr_r = A[3]
        self.start_r = int(A[4])
        self.end_r = int(A[5])

        self.sv_type = A[6]

        self.evidence_size = A[7]

        self.p_l = [float(x) for x in A[8].split(',')]
        self.p_r = [float(x) for x in A[9].split(',')]

    def __str__(self):
        return '\t'.join([self.chr_l, \
                         str(self.start_l), \
                         str(self.end_l), \
                         self.chr_r, \
                         str(self.start_r), \
                         str(self.end_r), \
                         self.sv_type, \
                         str(self.evidence_size)])

# I has 3 components [[start],[end],[p array]]
def align_intervals(I):
    start = -1
    end = -1
    new_I = []

    START = 0
    END = 1
    P = 2

    # find ends
    for i in I:
        if start == -1:
            start = i[START]
            end = i[END]
        else:
            if i[START] < start:
                start = i[START]

            if i[END] > end:
                end = i[END]

    for i in I:
        new_i = i[P]

        if i[START] > start:
            n = i[START] - start
            new_i = [0]*n + new_i

        if i[END] < end:
            n = end - i[END] 
            new_i = new_i + [0]*n
        
        new_I.append(new_i)
        
    return [start, end, new_I]
            

def cluster(BP):
    R = []
    L = []
    for b in BP:
        L.append([b.start_l,b.end_l,b.p_l])
        R.append([b.start_r,b.end_r,b.p_r])

    [start_R, end_R, a_R] = align_intervals(R)
    [start_L, end_L, a_L] = align_intervals(L)

    
    sum_L = [sum(x) for x in zip(*a_L)]
    sum_R = [sum(x) for x in zip(*a_R)]


    max_i_L = sum_L.index(max(sum_L))
    max_i_R = sum_R.index(max(sum_R))

    max_o_L = max_i_L + start_L
    max_o_R = max_i_R + start_R


    in_BP = []
    out_BP = []

    for b_i in range(len(BP)):
        if max_o_L >= BP[b_i].start_l and \
           max_o_L <= BP[b_i].end_l and \
           max_o_R >= BP[b_i].start_r and \
           max_o_R <= BP[b_i].end_r:
            in_BP.append(b_i)
        else:
            out_BP.append(b_i)

    for b_i in in_BP:


parser = OptionParser()

parser.add_option("-d",
    "--data_file",
    dest="data_file",
    help="Sorted lumpy data file")

(options, args) = parser.parse_args()

if not options.data_file:
    parser.error('Data not given')

f = open(options.data_file,'r')

BP_l = []
BP_starts_l = []
BP_ends_l = []

for l in f:
    bp = breakpoint(l)

    if (len(BP_l) == 0) or \
        (sum(BP_starts_l)/len(BP_starts_l) <= bp.end_l and \
         sum(BP_ends_l)/len(BP_ends_l) >= bp.start_l):
        BP_l.append(bp)
        BP_starts_l.append(bp.start_l)
        BP_ends_l.append(bp.end_l)
    else:
        # at this point all of the left ends intersect, but there could be
        # multiple clusters here that have different left ends so reprocess 
        # BP_l and cluster on the right side
        BP_r = []
        BP_starts_r = []
        BP_ends_r = []
        for bp_l in BP_l:
            if (len(BP_r) == 0) or \
                (sum(BP_starts_r)/len(BP_starts_r) <= bp_l.end_l and \
                 sum(BP_ends_r)/len(BP_ends_r) >= bp_l.start_l):
                BP_r.append(bp_l)
                BP_starts_r.append(bp_l.start_l)
                BP_ends_r.append(bp_l.end_l)
            else:
                cluster(BP_r)
                BP_r = []
                BP_starts_r = []
                BP_ends_r = []
                BP.append(bp_l)
                BP_starts_r.append(bp_l.start_l)
                BP_ends_r.append(bp_l.end_l)
        if len(BP_r) > 0:
            cluster(BP_r)

        BP = []
        BP_starts_l = []
        BP_ends_l = []
        BP.append(bp)
        BP_starts_l.append(bp.start_l)
        BP_ends_l.append(bp.end_l)

if len(BP) > 0:
    BP_r = []
    BP_starts_r = []
    BP_ends_r = []
    for bp_l in BP_l:
        if (len(BP_r) == 0) or \
            (sum(BP_starts_r)/len(BP_starts_r) <= bp_l.end_l and \
             sum(BP_ends_r)/len(BP_ends_r) >= bp_l.start_l):
            BP_r.append(bp_l)
            BP_starts_r.append(bp_l.start_l)
            BP_ends_r.append(bp_l.end_l)
        else:
            cluster(BP_r)
            BP_r = []
            BP_starts_r = []
            BP_ends_r = []
            BP.append(bp_l)
            BP_starts_r.append(bp_l.start_l)
            BP_ends_r.append(bp_l.end_l)
    if len(BP_r) > 0:
        cluster(BP_r)


