import sys
from sets import Set
import re

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

def parse_vcf(vcf_file_name, vcf_lines, vcf_headers, add_sname=True):
    header = ''
    samples = ''

    f = open(vcf_file_name, 'r')

    for l in f:
        if l[0] == '#':
            if l[1] != '#':
                samples = l.rstrip().split('\t')[9:]
            else:
                # ignore fileDate
                if l[:10] == '##fileDate':
                    continue
                if l not in vcf_headers:
                    vcf_headers.append(l)
        else:
            A = l.split('\t')
            if not 'SECONDARY' in A[7]:

                if add_sname and (samples != ''):
                    A[7] += ';' + 'SNAME=' + ','.join(samples)
                    l = '\t'.join(A)

                
                if 'SVTYPE=BND' in A[7]:
                    m = re.search(r"(\[|\])(.*)(\[|\])",A[4])
                    o_chr,o_pos = m.group(2).split(':')

                    if (o_chr == A[0]) and (('--:' in A[7]) != ('++' in A[7])):
                        neg_s = A[7].find('--:')
                        pos_s = A[7].find('++:')

                        if neg_s > 0:
                            neg_e = neg_s + A[7][neg_s:].find(';') 
                            pre=A[7][:neg_s]
                            mid=A[7][neg_s:neg_e]
                            post=A[7][neg_e:]
                            A[7] = pre + '++:0,' + mid + post
                        else:
                            pos_e = pos_s + A[7][pos_s:].find(';') 
                            pre=A[7][:pos_s]
                            mid=A[7][pos_s:pos_e]
                            post=A[7][pos_e:]
                            A[7] = pre + mid + ',--:0' + post

                        A[7] = 'SVTYPE=INV' + A[7][10:] + ';END=' + o_pos
                        A[4] = '<INV>'
                        vcf_lines.append('\t'.join(A))
                    else:
                        vcf_lines.append(l)
                else:
                    vcf_lines.append(l)

    return samples

def split_v(l):
    A = l.split('\t')
    m = to_map(A[7])

    chr_l = A[0]
    pos_l = int(A[1])

    chr_r = A[0]
    pos_r = int(A[1])
    if m['SVTYPE'] == 'BND':
        sep = '['
        if not sep in A[4]:
            sep = ']'
        s,e = [x for x in find_all(A[4],sep)]
        chr_r,pos_r = A[4][s+1:e].split(':')
        m['END'] = pos_r
        pos_r = int(pos_r)
    else:
        pos_r = int(m['END'])

    start_l = pos_l + int(m['CIPOS'].split(',')[0])
    end_l = pos_l + int(m['CIPOS'].split(',')[1])

    start_r = pos_r + int(m['CIEND'].split(',')[0])
    end_r = pos_r + int(m['CIEND'].split(',')[1])
        
    strands = m['STRANDS']

    return [m['SVTYPE'],chr_l,chr_r,strands,start_l,end_l,start_r,end_r,m]

def to_map(s):
    m = {}
    for k_v in s.split(';'):
        A = k_v.split('=')
        if len(A) > 1:
            m[A[0]] = A[1]
        else:
            m[A[0]] = None

    return m

def vcf_line_cmp(l1, l2):
    v1 = split_v(l1)
    v2 = split_v(l2)

    v1[3] = v1[3][:2]
    v2[3] = v2[3][:2]

    for i in range(len(v1)-1):
        if v1[i] != v2[i]:
            return cmp(v1[i],v2[i])
    return 0

def header_line_cmp(l1, l2):
    order = ['##source', \
             '##INFO', \
             '##ALT', \
             '##FORMAT',\
             '##SAMPLE']

    # make sure ##fileformat is first
    if l1[:12] == '##fileformat':
        return -1

    if l2[:12] == '##fileformat':
        return 1 

    # make sure #CHROM ... is last
    if l1[1] != '#':
        return 1
    elif l2[1] != '#':
        return -1

    if l1.find('=') == -1:
        return -1 
    if l2.find('=') == -1:
        return 1

    h1 = l1[:l1.find('=')]
    h2 = l2[:l2.find('=')]
    if h1 not in order:
        return -1 
    if h2 not in order:
        return 1
    return cmp(order.index(h1),order.index(h2))

class breakpoint:
    chr_l = ''
    start_l = 0
    end_l = 0
    p_l = []

    chr_r = ''
    start_r = 0
    end_r = 0
    p_r = []

    sv_type = ''

    strands = ''
    
    l = ''

    def __init__(self, 
                 l,
                 percent_slop=0,
                 fixed_slop=0):
        self.l = l

        [self.sv_type,\
        self.chr_l, \
        self.chr_r,\
        self.strands,
        self.start_l,\
        self.end_l,\
        self.start_r, \
        self.end_r, 
        m] = split_v(l)

        self.p_l = [float(x) for x in m['PRPOS'].split(',')]
        self.p_r = [float(x) for x in m['PREND'].split(',')]

        slop_prob = 1e-100
        if ((percent_slop > 0) or (fixed_slop > 0)):

            l_slop = int(max(percent_slop*(self.end_l-self.start_l),fixed_slop))
            r_slop = int(max(percent_slop*(self.end_r-self.start_r),fixed_slop))

            # pad each interval with slop_prob on each side.
            self.start_l = self.start_l-l_slop
            self.end_l = self.end_l+l_slop
            new_p_l = [slop_prob] * l_slop + self.p_l + [slop_prob] * l_slop

            self.start_r = self.start_r-r_slop
            self.end_r = self.end_r+r_slop
            new_p_r = [slop_prob] * r_slop + self.p_r + [slop_prob] * r_slop

            # chew off overhang if self.start_l or self.start_r less than 0
            if self.start_l < 0:
                new_p_l = new_p_l[-self.start_l:]
                self.start_l = 0
            if self.start_r < 0:
                new_p_r = new_p_r[-self.start_r:]
                self.start_r = 0

            # normalize so each probability curve sums to 1.
            sum_p_l = sum(new_p_l)
            self.p_l = [float(x)/sum_p_l for x in new_p_l]
            sum_p_r = sum(new_p_r)
            self.p_r = [float(x)/sum_p_r for x in new_p_r]

            # old_l = float(self.end_l - self.start_l + 1)
            
            # self.start_l = max(0,self.start_l-l_slop)
            # self.end_l = self.end_l+l_slop

            # new_l = float(self.end_l - self.start_l + 1)

            # new_p_l = []
            # for i in range(self.end_l-self.start_l+1):
            #     p = i/new_l
            #     old_i = int(p*old_l)
            #     new_p_l.append(self.p_l[old_i])
            # sum_p_l = sum(new_p_l)
            # self.p_l = [float(x)/sum_p_l for x in new_p_l]

            # old_r = float(self.end_r - self.start_r + 1)

            # self.start_r = max(0,self.start_r-r_slop)
            # self.end_r = self.end_r+r_slop

            # new_r = float(self.end_r - self.start_r + 1)

            # new_p_r = []
            # for i in range(self.end_r-self.start_r+1):
            #     p = float(i)/new_r
            #     old_i = int(p*old_r)
            #     new_p_r.append(self.p_r[old_i])
            # sum_p_r = max(1,sum(new_p_r))
            # self.p_r = [float(x)/sum_p_r for x in new_p_r]

    def __str__(self):
        return '\t'.join([str(x) for x in [self.chr_l, \
                                           self.start_l,\
                                           self.end_l,\
                                           self.chr_r,\
                                           self.start_r, \
                                           self.end_r, 
                                           self.sv_type,\
                                           self.strands,\
                                           self.p_l,
                                           self.p_r]])
    def ovl(self, b):
        if (self.chr_l != b.chr_l) or \
            (self.chr_r != b.chr_r) or \
            (self.sv_type != b.sv_type):
            return 0
        #get left common interval
        c_start_l, c_end_l = [max(self.start_l, b.start_l), \
                              min(self.end_l, b.end_l)]
        #get right common interval
        c_start_r, c_end_r = [max(self.start_r, b.start_r), \
                              min(self.end_r, b.end_r)]

        c_l_len = c_end_l - c_start_l + 1
        c_r_len = c_end_r - c_start_r + 1

        if (c_l_len < 1) or (c_r_len < 1):
            return 0

        self_start_off_l = c_start_l - self.start_l
        b_start_off_l = c_start_l - b.start_l

        self_start_off_r = c_start_r - self.start_r
        b_start_off_r = c_start_r - b.start_r

        ovl_l = 0
        for i in range(c_l_len):
            ovl_l += min(self.p_l[i+self_start_off_l], b.p_l[i+b_start_off_l])

        ovl_r = 0
        for i in range(c_r_len):
            ovl_r += min(self.p_r[i+self_start_off_r], b.p_r[i+b_start_off_r])

        return ovl_l * ovl_r

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


class node:
    b = None
    color = -1
    edges = None
    def __init__(self, b):
        self.b = b
        self.color = -1
        self.edges = []

def connect(G, B, t):
    #first we need to add all of the elements in B to the graph

    # each node in the graph has 3 elements
    # 0: breakpoint
    # 1: color
    # 2: list of edges
    #    each edge has 2 elements
    #    0: correspondined node id in G
    #    1: weight (ovl score)

    b_ids = []

    for b in B:
        next_id = len(G)
        b_ids.append(next_id)
        #G[next_id] = [b, -1, []]
        G[next_id] = node(b)

    for i in range(len(B)):
        for j in range(len(B)):
            if i != j:
                ovl = B[i].ovl(B[j])
                if ovl > t:
                    #G[b_ids[i]][2].append([b_ids[j], ovl])
                    G[b_ids[i]].edges.append([b_ids[j], ovl])

def bron_kerbosch(G, R, P, X):
    if (len(P) == 0) and (len(X) == 0):
        yield R
    for v in P:
        V = Set([v])
        N = Set([g[0] for g in G[v].edges])
    
        for r in bron_kerbosch(G, \
                               R.union(V), \
                               P.intersection(N), 
                               X.intersection(N)):
            yield r

        P = P - V
        X = X.union(V)




