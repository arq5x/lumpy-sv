#!/usr/bin/env python -u
from operator import add
import time
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

    v_id = 0
    score = 0
    ids_str = ''
    strands_str = ''
    max_str = ''
    nine_five_str = ''
    sv_type = ''
    evidence_size = 0
    src_file = ''

    def __init__(self, l):
        A = l.rstrip().split('\t')
        self.chr_l = A[0]
        self.start_l = int(A[1])
        self.end_l = int(A[2])

        self.chr_r = A[3]
        self.start_r = int(A[4])
        self.end_r = int(A[5])

        self.v_id = A[6]
        self.score = A[7]
        self.strand_l = A[8]
        self.strand_r = A[9]

        self.sv_type = A[10]

        #for ec in A[11][4:].split(';'):
            #e,c = ec.split(',')
            #self.evidence_size += int(c)

        self.ids_str = A[11]
        self.strands_str = A[12]
        self.max_str = A[13]
        self.nine_five_str = A[14]

        self.p_l = [float(x) for x in A[15][3:].split(',')]
        self.p_r = [float(x) for x in A[16][3:].split(',')]

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

def trim(A):
    clip_start = 0
    for i in range(len(A)):
        if A[i] == 0:
            clip_start += 1
        else:
            break

    clip_end = 0
    for i in range(len(A)-1,-1,-1):
        if A[i] == 0:
            clip_end += 1
        else:
            break               
    return [clip_start, clip_end]

def max_filter(BP):
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

    in_BP_i = []
    out_BP_i = []

    for b_i in range(len(BP)):
        if max_o_L >= BP[b_i].start_l and \
           max_o_L <= BP[b_i].end_l and \
           max_o_R >= BP[b_i].start_r and \
           max_o_R <= BP[b_i].end_r:
            in_BP_i.append(b_i)
        else:
            out_BP_i.append(b_i)

    return [in_BP_i, out_BP_i]

def merge(BP):

    # optionally filter the collection
    #in_BP_i, out_BP_i = max_filter(BP)

    # no filter
    in_BP_i = range(len(BP))
    out_BP_i = []

    if len(in_BP) > 0:
        # get the common interval
        max_start_l = 0
        min_end_l =  sys.maxint

        max_start_r = 0
        min_end_r =  sys.maxint

        for b_i in in_BP:
            # left side
            if BP[i].start_l > max_start_l:
                max_start_l = BP[i].start_l
            if BP[i].end_l < min_end_l:
                min_end_l = BP[i].end_l

            # right side
            if BP[i].start_r > max_start_r:
                max_start_r = BP[i].start_r

            if BP[i].end_r < min_end_r:
                min_end_r = BP[i].end_r



    if len(out_BP) > 0:
        re_merge = []
        for i in out_BP:
            re_merge.append(BP[i])

        if len(in_BP) == 0:
            print
            print '\t', max_o_L, max_o_R
            for r in re_merge:
                print r
            exit(1)
        else:
            print 'R\t',
            re_merge.sort(key=lambda x: x.start_l)
            re_merge.sort(key=lambda x: x.sv_type)
            re_merge.sort(key=lambda x: x.chr_l)
            r_cluster(re_merge)


#        print '\t'.join([\
#                chr_L, \
#                str(new_start_L), \
#                str(new_end_L), \
#                chr_R, \
#                str(new_start_R), \
#                str(new_end_R), \
#                "0", \
#                str(len(in_BP)), \
#                BP[0].strand_l, \
#                BP[0].strand_r, \
#                BP[0].sv_type, \
#                "LP:" + ','.join(\
#                    [str(x) for x in p_L[clip_start_L:len(p_L)-clip_end_L]]), \
#                "RP:" + ','.join(\
#                    [str(x) for x in p_R[clip_start_R:len(p_R)-clip_end_R]])\
#                ])

# Scan the list of calls that intersect on the left side.  These must first be
# esolted based on the rigth side coordinates.  Merge when a call does not
# intersect the current collection.
def r_cluster(BP_l):
    # need to resort based on the right side, then extract clusters
    BP_l.sort(key=lambda x: x.start_r)
    BP_l.sort(key=lambda x: x.sv_type)
    BP_l.sort(key=lambda x: x.chr_r)

    BP_r = []
    BP_starts_r = []
    BP_ends_r = []
    for bp_l in BP_l:
        if (len(BP_r) == 0) or \
            (sum(BP_starts_r)/len(BP_starts_r) <= bp_l.end_r and \
             sum(BP_ends_r)/len(BP_ends_r) >= bp_l.start_r):
            BP_r.append(bp_l)
            BP_starts_r.append(bp_l.start_r)
            BP_ends_r.append(bp_l.end_r)
        else:    
            merge(BP_r)
            BP_r = []
            BP_starts_r = []
            BP_ends_r = []
            BP_r.append(bp_l)
            BP_starts_r.append(bp_l.start_r)
            BP_ends_r.append(bp_l.end_r)
    if len(BP_r) > 0:
        merge(BP_r)

# The input file should be sorted, collect calls that intersect on the left
# side as soon as a call is found that does not intersect the current
# collection repeat the process for the right side
def l_cluster(in_file):

    f = open(in_file, 'r')

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

            r_cluster(BP_l)

            BP_l = []
            BP_starts_l = []
            BP_ends_l = []
            BP_l.append(bp)
            BP_starts_l.append(bp.start_l)
            BP_ends_l.append(bp.end_l)

    if len(BP_l) > 0:
        r_cluster(BP_l)

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():

    usage = """%prog -i <file>
lmerge
Author: Ryan Layer & Ira Hall
Description: merge lumpy calls.
Version: ira_7
    """

    parser = OptionParser(usage)
    parser.add_option("-i", \
                      "--inFile", \
                      dest="inFile",
                      help="A sorted lumpy output file generated by " + \
                           "lsort; or stdin (-i stdin). Column 7 must " + \
                           "have the format sample:variantID", \
                           metavar="FILE")

#    parser.add_option("-c", \
#                      "--configFile", \
#                      dest="configFile", \
#                      help="The config file: sample,id,evidenceType,path", \
#                      metavar="FILE")

    (opts, args) = parser.parse_args()

    #if opts.inFile is None or opts.configFile is None:
    if opts.inFile is None:
        parser.print_help()
        print
    else:
        try:
            #l_cluster(opts.inFile, opts.configFile)
            l_cluster(opts.inFile)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return

if __name__ == "__main__":
    sys.exit(main())
