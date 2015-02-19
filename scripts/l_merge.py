#!/usr/bin/env python -u
from operator import add
import time
import sys
import numpy as np
import glob
from operator import add
from optparse import OptionParser
from sets import Set
import l_bp

def merge(BP):

    if len(BP) == 1:
        ##tack on id to SNAME
        A = BP[0].l.rstrip().split('\t')
        A[7]+= ':' + A[2]
        print '\t'.join(A)
        return 0

    G = {}
    l_bp.connect(G, BP, 0)
    #print G[0].b.chr_l
    #for i in G:
        #print i,[g[0] for g in G[i].edges]

    C = []
    _G = G.copy()
    while len(_G) != 0:
        R = Set()
        X = Set()
        P = Set(_G.keys())
        clique = [x for x in l_bp.bron_kerbosch(_G, R, P, X)]
        max_clique = sorted(clique, key=len)[0]
        C.append(list(max_clique))
        # remove these from the graph
        for g in _G:
            E = [e for e in _G[g].edges if e[0] not in clique]
            G[g].edges = E
        for c in max_clique:
            del _G[c]
    
    for c in C:
        #print len(c)
        L = []
        R = []
        for g_i in c:
            b = G[g_i].b
            L.append([b.start_l,b.end_l,b.p_l])
            R.append([b.start_r,b.end_r,b.p_r])

        [start_R, end_R, a_R] = l_bp.align_intervals(R)
        [start_L, end_L, a_L] = l_bp.align_intervals(L)

        p_L = [1] * len(a_L[0])
        p_R = [1] * len(a_R[0])

        for c_i in range(len(c)):
            for i in range(len(a_L[c_i])):
                p_L[i] = p_L[i] * a_L[c_i][i]

            for i in range(len(a_R[c_i])):
                p_R[i] = p_R[i] * a_R[c_i][i]


        [clip_start_L, clip_end_L] = l_bp.trim(p_L)
        [clip_start_R, clip_end_R] = l_bp.trim(p_R)

        new_start_L = start_L + clip_start_L
        new_end_L = end_L - clip_end_L

        new_start_R = start_R + clip_start_R
        new_end_R = end_R - clip_end_R

        p_L = p_L[clip_start_L:len(p_L)-clip_end_L]
        p_R = p_R[clip_start_R:len(p_R)-clip_end_R]

        s_p_L = sum(p_L)
        s_p_R = sum(p_R)

        p_L = [x/s_p_L for x in p_L]
        p_R = [x/s_p_R for x in p_R]

        max_i_L = p_L.index(max(p_L))
        max_i_R = p_R.index(max(p_R))

        ninefive_i_L_start = max_i_L
        ninefive_i_L_end = max_i_L
        ninefive_i_L_total = p_L[max_i_L]
        updated = 0
        while (ninefive_i_L_total < 0.95):
                if ninefive_i_L_start > 0:
                    ninefive_i_L_start -= 1
                    ninefive_i_L_total += p_L[ninefive_i_L_start]
                    updated = 1
                if ninefive_i_L_end < len(p_L)-1:
                    ninefive_i_L_end += 1
                    ninefive_i_L_total += p_L[ninefive_i_L_end]
                    updated = 1
                if not updated:
                    break

        ninefive_i_R_start = max_i_R
        ninefive_i_R_end = max_i_R
        ninefive_i_R_total = p_R[max_i_R]
        updated = 0
        while (ninefive_i_R_total < 0.95):
                if ninefive_i_R_start > 0:
                    ninefive_i_R_start -= 1
                    ninefive_i_R_total += p_R[ninefive_i_R_start]
                    updated = 1
                if ninefive_i_R_end < len(p_R)-1:
                    ninefive_i_R_end += 1
                    ninefive_i_R_total += p_R[ninefive_i_R_end]
                    updated = 1
                if not updated:
                    break
 
        CIPOS95=str(ninefive_i_L_start) + ',' + str(ninefive_i_L_end)
        CIEND95=str(ninefive_i_R_start) + ',' + str(ninefive_i_R_end)

        CHROM = G[c[0]].b.chr_l
        POS = new_start_L + max_i_L
        ID = 1
        REF = 'N'

        ALT = ''
        if G[c[0]].b.sv_type == 'BND':
            ALT = 'N]' + G[c[0]].b.chr_r + ':' + str(new_start_R + max_i_R)
        else:
            ALT = '<' + G[c[0]].b.sv_type + '>'
        QUAL = '.'
        FILTER = '.'
        FORMAT = G[c[0]].b.l.split('\t')[8]
        SVTYPE = G[c[0]].b.sv_type

        STRANDS = ''
        strand_map = {}
        e_type_map = {}

        SU = 0
        PE = 0
        SR = 0

        s_name_list = []

        gt_list = [] 

        for g_i in c:
            A = G[g_i].b.l.rstrip().split('\t')
            m = l_bp.to_map(A[7])


            for strand_entry in m['STRANDS'].split(','):
                s_type,s_count = strand_entry.split(':')
                if s_type not in strand_map:
                    strand_map[s_type] = 0
                strand_map[s_type] += int(s_count)

            SU += int(m['SU'])
            PE += int(m['PE'])
            SR += int(m['SR'])

            s_name_list.append(m['SNAME'] + ':' + A[2])

            gt_list += A[9:]

        SNAME=','.join(s_name_list)

        GTS = '\t'.join(gt_list)

        strand_types_counts = []
        for strand in strand_map:
            strand_types_counts.append(strand + ':' + str(strand_map[strand]))
        STRANDS = ','.join(strand_types_counts)

        SVLEN = (new_start_L + max_i_L) - (new_start_R + max_i_R)
        END = new_start_R + max_i_R
        CIPOS=','.join([str(x) for x in [-1*max_i_L, len(p_L) - max_i_L - 1]])
        CIEND=','.join([str(x) for x in [-1*max_i_R, len(p_R) - max_i_R - 1]])
        IMPRECISE='IMPRECISE'
        PRPOS=','.join([str(x) for x in p_L])
        PREND=','.join([str(x) for x in p_R])

        INFO = ';'.join(['SVTYPE='   + str(SVTYPE),
                         'STRANDS='  + str(STRANDS),
                         'END='      + str(END),
                         'CIPOS='    + str(CIPOS),
                         'CIEND='    + str(CIEND),
                         'CIPOS95='  + str(CIPOS95),
                         'CIEND95='  + str(CIEND95),
                                       str(IMPRECISE),
                         'SU='       + str(SU),
                         'PE='       + str(PE),
                         'SR='       + str(SR),
                         'PRPOS='    + str(PRPOS),
                         'PREND='    + str(PREND),
                         'SNAME='    + str(SNAME)])

        O = [CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GTS]

        print '\t'.join([str(o) for o in O])

def r_cluster(BP_l):
    # need to resort based on the right side, then extract clusters
    BP_l.sort(key=lambda x: x.start_r)
    BP_l.sort(key=lambda x: x.chr_r)

    BP_r = []
    BP_max_end_r = -1
    BP_chr_r = ''

    for b in BP_l:
        if (len(BP_r) == 0) or \
           ((b.start_r <= BP_max_end_r) and \
           (b.chr_r == BP_chr_r)):
            BP_r.append(b)
            BP_max_end_r = max(BP_max_end_r, b.end_r)
            BP_chr_r = b.chr_r
        else:
            #print len(BP_r)
            merge(BP_r)
            BP_r = [b]
            BP_max_end_r = b.end_r
            BP_chr_r = b.chr_r
 
    if len(BP_r) > 0:
        #print len(BP_r)
        merge(BP_r)

def l_cluster(file_name):
    vcf_lines = []
    vcf_headers = Set()
    r = l_bp.parse_vcf(file_name, vcf_lines, vcf_headers)

    BP_l = []
    BP_sv_type = ''
    BP_max_end_l = -1
    BP_chr_l = ''

    for l in vcf_lines:
        b = l_bp.breakpoint(l)

        if (len(BP_l) == 0) or \
           ((b.start_l <= BP_max_end_l) and \
            (b.chr_l == BP_chr_l) and \
            (b.sv_type == BP_sv_type)):
            BP_l.append(b)
            BP_max_end_l = max(BP_max_end_l, b.end_l)
            BP_chr_l = b.chr_l
            BP_sv_type = b.sv_type
        else:
            #print len(BP_l)
            r_cluster(BP_l)
            BP_l = [b]
            BP_max_end_l = b.end_l
            BP_sv_type = b.sv_type
            BP_chr_l = b.chr_l

    if len(BP_l) > 0:
        #print len(BP_l)
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

    (opts, args) = parser.parse_args()

    #if opts.inFile is None or opts.configFile is None:
    if opts.inFile is None:
        parser.print_help()
        print
    else:
        try:
            l_cluster(opts.inFile)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return

if __name__ == "__main__":
    sys.exit(main())
