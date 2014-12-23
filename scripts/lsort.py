#!/usr/bin/env python
import sys
import numpy as np
import glob
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

    src_file = ''


    def __init__(self, line, src_file):
        self.src_file = src_file

        A = l.rstrip().split('\t')
        self.chr_l = A[0]
        self.start_l = int(A[1])
        self.end_l = int(A[2])

        self.chr_r = A[3]
        self.start_r = int(A[4])
        self.end_r = int(A[5])

        self.sv_type = A[10][5:]

        for ec in A[11][4:].split(';'):
            e,c = ec.split(',')
            self.evidence_size += int(c)

        #self.p_l = [float(x) for x in A[15][3:].split(',')]
        #self.p_r = [float(x) for x in A[16][3:].split(',')]
        self.p_l = A[15][3:]
        self.p_r = A[16][3:]

    def __str__(self):
        return '\t'.join([self.chr_l, \
                         str(self.start_l), \
                         str(self.end_l), \
                         self.chr_r, \
                         str(self.start_r), \
                         str(self.end_r), \
                         self.sv_type, \
                         str(self.evidence_size),
                         self.p_l,\
                         self.p_r ]) 
    def __cmp__(self,other):
       
        presidence = [ 
                       [self.sv_type,other.sv_type],\
                       [self.chr_l,other.chr_l], \
                       [self.chr_r,other.chr_r], \
                       [self.start_l,other.start_l], \
                       [self.end_l,other.end_l], \
                       [self.start_r,other.end_l], 
                       ]        

        for p in presidence:
            if p[0] != p[1]:
                return cmp(p[0],p[1])

        return 0


parser = OptionParser()

parser.add_option("-d",
    "--data_file_dir",
    dest="data_file_dir",
    help="Data file directory")

(options, args) = parser.parse_args()

if not options.data_file_dir:
    parser.error('Data file directory not given')

BP = []

for file_name in glob.glob(options.data_file_dir):
    f = open(file_name,'r')

    for l in f:
        BP.append(breakpoint(l,file_name))
    f.close()

BP.sort()

for bp in BP:
    print bp
