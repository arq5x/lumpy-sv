#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>

#include "htslib/sam.h"
#include "htslib/bgzf.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#define MIN_NON_OVERLAP 20
#define MAX_UNMAPPED_BASES 50
#define MIN_INDEL_SIZE 50

struct line {
    uint32_t pos, raLen, rapos, qaLen, sclip, eclip, SQO, EQO;
    char strand;
    char *chrm;
};


// from Dave Larson's extract_sv_reads
static inline size_t _pre_seq_bytes(bam1_t const* b) {
	return (b->core.n_cigar<<2) + b->core.l_qname;
}
static inline size_t _seq_qual_bytes(bam1_t const* b) {
	// sequence bytes + quality bytes
	return (((b->core.l_qseq + 1)>>1) + b->core.l_qseq);
}
static inline size_t _post_qual_bytes(bam1_t const* b) {
	// Total - pre-sequence bytes - sequence bytes - quality bytes
	return (b->l_data - _pre_seq_bytes(b) - _seq_qual_bytes(b));
}
bam1_t *drop_seq_qual(bam1_t *_aln) {
	bam1_t *b = bam_dup1(_aln);
	uint8_t* data = _aln->data;
	int m_data = _aln->m_data;
	// NOTE Copying everything BUT sequence and quality
	memcpy(data, b->data, _pre_seq_bytes(b));
	memcpy(data + _pre_seq_bytes(b), bam_get_aux(b), _post_qual_bytes(b));
	*_aln = *b;
	_aln->l_data = b->l_data - _seq_qual_bytes(b);
	_aln->core.l_qseq = 0;
	_aln->m_data = m_data;
	_aln->data = data;
	bam_destroy1(b);
	return _aln;
}

// This will parse a base 10 int, and change ptr to one char beyond the end of the number.
int parseNextInt(char **ptr)
{
    int num = 0;
    char curChar;
    for (curChar = (*ptr)[0]; curChar != 0; curChar = (++(*ptr))[0])
    {
        int digit = curChar - '0';
        if (digit >= 0 && digit <= 9) num = num*10 + digit;
        else break;
    }
    return num;
}

// This will the current char, and move the ptr ahead by one.
char parseNextOpCode(char **ptr)
{
    return ((*ptr)++)[0];
}

// This just test for end of string.
bool moreCigarOps(char *ptr)
{
    return (ptr[0] != 0);
}

void calcAlnOffsets(uint32_t *cigar,
                    uint32_t n_cigar,
                    uint32_t sa_pos,
                    char sa_strand,
                    struct line *l)
{
    l->raLen = 0;
    l->qaLen = 0;
    l->sclip = 0;
    l->eclip = 0;
    l->SQO = 0;
    l->EQO = 0;
    bool first = true;

    uint32_t k;
    for (k = 0; k < n_cigar; ++k)
    {
        uint32_t opLen = bam_cigar_oplen(cigar[k]);
        char opCode = bam_cigar_opchr(cigar[k]);
        if      (opCode == 'M' || opCode == '=' || opCode == 'X')
        {
            l->raLen += opLen;
            l->qaLen += opLen;
            first = false;
        }
        else if (opCode == 'S' || opCode == 'H')
        {
            if (first) l->sclip += opLen;
            else       l->eclip += opLen;
        }
        else if (opCode == 'D' || opCode == 'N')
        {
            l->raLen += opLen;
        }
        else if (opCode == 'I')
        {
            l->qaLen += opLen;
        }
    }
    //*rapos = str2pos(line->fields[POS]);
    l->rapos = sa_pos;
    if (sa_strand == '+')
    {
        l->pos = l->rapos - l->sclip;
        l->SQO = l->sclip;
        l->EQO = l->sclip + l->qaLen - 1;
    }
    else
    {
        l->pos = l->rapos + l->raLen + l->eclip - 1;
        l->SQO = l->eclip;
        l->EQO = l->eclip + l->qaLen - 1;
    }
}



void calcOffsets(char *cigar,
                 uint32_t sa_pos,
                 char sa_strand,
                 struct line *l)
{
    l->raLen = 0;
    l->qaLen = 0;
    l->sclip = 0;
    l->eclip = 0;
    l->SQO = 0;
    l->EQO = 0;
    bool first = true;
    while (moreCigarOps(cigar))
    {
        int opLen = parseNextInt(&cigar);
        char opCode = parseNextOpCode(&cigar);
        if      (opCode == 'M' || opCode == '=' || opCode == 'X')
        {
            l->raLen += opLen;
            l->qaLen += opLen;
            first = false;
        }
        else if (opCode == 'S' || opCode == 'H')
        {
            if (first) l->sclip += opLen;
            else       l->eclip += opLen;
        }
        else if (opCode == 'D' || opCode == 'N')
        {
            l->raLen += opLen;
        }
        else if (opCode == 'I')
        {
            l->qaLen += opLen;
        }
    }
    //*rapos = str2pos(line->fields[POS]);
    l->rapos = sa_pos;
    if (sa_strand == '+')
    {
        l->pos = l->rapos - l->sclip;
        l->SQO = l->sclip;
        l->EQO = l->sclip + l->qaLen - 1;
    }
    else
    {
        l->pos = l->rapos + l->raLen + l->eclip - 1;
        l->SQO = l->eclip;
        l->EQO = l->eclip + l->qaLen - 1;
    }
}


void split_sa_tag(char *sa_tag,
                  char **chrm,
                  uint32_t *pos,
                  char *strand,
                  char **cigar)
{
    *chrm = strtok(sa_tag, ",");
    *pos = atoi(strtok(NULL, ","));
    *strand = strtok(NULL, ",")[0];
    *cigar = strtok(NULL, ",");
}

int count_tags(char *sa_tag)
{
    uint32_t i, c = 0;
    for (i = 0; i < strlen(sa_tag); ++i) {
        if (sa_tag[i] == ';')
            c += 1;
    }

    return c;
}

int main(int argc, char **argv)
{
	// without these 2 lines, htslib sometimes tries to download a part of the sequence
	// even though the -f reference was provided.
	setenv("REF_CACHE", "", 0);
	setenv("REF_PATH", "fake_value_so_no_download", 0);
    if (argc < 4)
        errx(1,
             "usage\t:%s -f <optional-reference> <bam> <split out> <discord out> (optional #threads)",
             argv[0]);


    int c;
    char *fasta = NULL;
    while (( c = getopt(argc, argv, "f:")) != -1) {
        fasta = optarg;
    }

    char *bam_file_name = argv[optind];
    char *split_file_name = argv[1+optind];
    char *disc_file_name = argv[2+optind];
    int threads = 2;
    if (argc == 4+optind) {
        threads = atoi(argv[3+optind]);
    }


    samFile *in = sam_open(bam_file_name, "rb");
    if(in == NULL) {
        errx(1, "Unable to open BAM/SAM file.");
    }
    // 0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x20 | 0x40 | 0x80 | 0x100 | 0x200 | 0x800
    // decode everything but base-qual
    hts_set_opt(in, CRAM_OPT_REQUIRED_FIELDS, 3071);
    hts_set_opt(in, CRAM_OPT_DECODE_MD, 1);
    if (fasta != NULL) {
        hts_set_fai_filename(in, fasta);
    }
    if (threads > 1) {
		hts_set_threads(in, threads);
    }

    bam_hdr_t *hdr = sam_hdr_read(in);

    samFile *disc = sam_open(disc_file_name, "wb1");
    samFile *split = sam_open(split_file_name, "wb1");

    int r = sam_hdr_write(disc, hdr);
    r = sam_hdr_write(split, hdr);

    bam1_t *aln = bam_init1();
    int ret;
	int dropped;

    while(ret = sam_read1(in, hdr, aln) >= 0) {
        if (((aln->core.flag) & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP)) != 0)
            continue;

		dropped = 0;
        if (((aln->core.flag) & (BAM_FPROPER_PAIR | BAM_FMUNMAP | BAM_FSUPPLEMENTARY)) == 0) {
			dropped = 1;
			aln = drop_seq_qual(aln);
            r = sam_write1(disc, hdr, aln);
		}

        uint8_t *sa = bam_aux_get(aln, "SA");
		if (sa == 0) {
			continue;
		}

        char *sa_tag = strdup(bam_aux2Z(sa));
        if ( count_tags(sa_tag) == 1) {
            char *chrm, strand, *cigar;
            uint32_t pos;
            split_sa_tag(sa_tag,
                         &chrm,
                         &pos,
                         &strand,
                         &cigar);

            struct line sa, al;

            calcOffsets(cigar,
                        pos,
                        strand,
                        &sa);
            sa.chrm = chrm;
            sa.strand = strand;


            calcAlnOffsets(bam_get_cigar(aln),
                           aln->core.n_cigar,
                           aln->core.pos + 1,
                           bam_is_rev(aln) ? '-' : '+',
                           &al);
            al.chrm = hdr->target_name[aln->core.tid];
            al.strand = bam_is_rev(aln) ? '-' : '+';

            struct line *left = &al, *right = &sa;

            if (left->SQO > right->SQO) {
                left = &sa;
                right = &al;
            }

            int overlap = MAX(1 + MIN(left->EQO, right->EQO) -
                    MAX(left->SQO, right->SQO), 0);
            int alen1 = 1 + left->EQO - left->SQO;
            int alen2 = 1 + right->EQO - right->SQO;
            int mno = MIN(alen1-overlap, alen2-overlap);
            if (mno < MIN_NON_OVERLAP)
                continue;

            if ( (strcmp(left->chrm, right->chrm) == 0) &&
                 (left->strand == right->strand) ) {

                int leftDiag, rightDiag, insSize;
                if (left->strand == '-') {
                    leftDiag = left->rapos - left->sclip;
                    rightDiag = (right->rapos + right->raLen) -
                            (right->sclip + right->qaLen);
                    insSize = rightDiag - leftDiag;
                } else {
                    leftDiag = (left->rapos + left->raLen) -
                            (left->sclip + left->qaLen);
                    rightDiag = right->rapos - right->sclip;
                    insSize = leftDiag - rightDiag;
                }
                int desert = right->SQO - left->EQO - 1;
                if ((abs(insSize) < MIN_INDEL_SIZE) ||
                    ((desert > 0) && (
                        (desert - (int)MAX(0, insSize)) >
                        MAX_UNMAPPED_BASES)))
                    continue;
            }

            char *qname =  bam_get_qname(aln);
            if ((aln->core.flag & 64) == 64)
                qname[0] = 'A';
            else
                qname[0] = 'B';

			if (dropped == 0) {
    			aln = drop_seq_qual(aln);
			}
            r = sam_write1(split, hdr, aln);
        }
        free(sa_tag);
    }
    bam_destroy1(aln);
    sam_close(disc);
    sam_close(split);
    bam_hdr_destroy(hdr);
    sam_close(in);
    if(ret < -1) {
        errx(1, "lumpy_filter: error reading bam: %s\n", bam_file_name);
    }
}
