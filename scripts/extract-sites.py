from __future__ import print_function
import sys
import gzip
import collections

def xopen(f):
    return (gzip.open if f.endswith(".gz") else open)(f)

sites = collections.OrderedDict()
files = sys.argv[1:]

def key(toks):
    return (toks[0], toks[1], toks[3], toks[4])

def anynonref(gts):
    return any(gt.startswith(("0/1", "1/1")) for gt in gts)


header = []
for f in files:
    header = []
    for toks in (l.rstrip().split("\t") for l in xopen(f)):
        if toks[0][0] == "#":
            if toks[0] == "#CHROM":
                toks = toks[:8]
            header.append("\t".join(toks))
            continue
        if anynonref(toks[9:]):
            sites[key(toks)] = "\t".join(toks[:8])

print("\n".join(header))
for k, line in sites.items():
    print(line)
