#!/usr/bin/env python
import sys
import numpy as np
import glob
from optparse import OptionParser
import l_bp
from sets import Set

def main():
    usage ="""%prog <VCF file 1> <VCF file 2> ... <VCF file N>

l_sort
Author: Ryan Layer, Colby Chiang, & Ira Hall
Description: sort N VCF files into a single file
Version: 0.01
"""

    if len(sys.argv) < 2:
        exit(1)

    vcf_file_names = sys.argv[1:]

    vcf_lines = []
    vcf_headers = Set()

    vcf_headers.add("##INFO=<ID=SVNAME,Number=.,Type=String," + \
            "Description=\"Source sample name\">\n")

    vcf_headers.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t" + \
            "VARIOUS\n")

    for vcf_file_name in vcf_file_names:
        r = l_bp.parse_vcf(vcf_file_name, vcf_lines, vcf_headers)

    for h in sorted(vcf_headers, cmp=l_bp.header_line_cmp):
        print h,

    vcf_lines.sort(cmp=l_bp.vcf_line_cmp)
    for v in vcf_lines:
        print v,

if __name__ == "__main__":
    sys.exit(main())
