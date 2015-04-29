#/*****************************************************************************
 * SV_Bedpe.cpp
 * (c) 2012 - Ryan M. Layer
 * Hall Laboratory
 * Quinlan Laboratory
 * Department of Computer Science
 * Department of Biochemistry and Molecular Genetics
 * Department of Public Health Sciences and Center for Public Health Genomics,
 * University of Virginia
 * rl6sf@virginia.edu
 *
 * Licenced under the GNU General Public License 2.0 license.
 * ***************************************************************************/
#include "BamAncillary.h"
using namespace BamTools;

#include "SV_BedpeReader.h"
#include "SV_BreakPoint.h"
#include "SV_Bedpe.h"
#include "log_space.h"
#include "bedFilePE.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <math.h>

using namespace std;

//{{{ SV_Bedpe:: SV_Bedpe(const BEDPE *bedpeEntry)
SV_Bedpe::
SV_Bedpe(const BEDPE *bedpeEntry,
         int _weight,
         int _ev_id,
         SV_BedpeReader *_reader)
{
    reader = _reader;

    ev_id = _ev_id;
    struct interval tmp_a, tmp_b;

    tmp_a.start = bedpeEntry->start1;
    tmp_a.end = bedpeEntry->end1 - 1;
    tmp_a.chr = bedpeEntry->chrom1;
    tmp_a.strand = bedpeEntry->strand1[0];

    tmp_b.start = bedpeEntry->start2;
    tmp_b.end = bedpeEntry->end2 - 1;
    tmp_b.chr = bedpeEntry->chrom2;
    tmp_b.strand = bedpeEntry->strand2[0];

    if ( tmp_a.chr.compare(tmp_b.chr) > 0 ) {
        side_l = tmp_a;
        side_r = tmp_b;
    } else if ( tmp_a.chr.compare(tmp_b.chr) < 0 ) {
        side_l = tmp_b;
        side_r = tmp_a;
    } else { // ==
        if (tmp_a.start > tmp_b.start) {
            side_l = tmp_b;
            side_r = tmp_a;
        } else {
            side_l = tmp_a;
            side_r = tmp_b;
        }
    }

    vector<string>::const_iterator it;
    type = -1;
    for (it = bedpeEntry->fields.begin();
            it != bedpeEntry->fields.end(); ++it) {

        if ( it->find("TYPE:") == 0 ) {
            string type_string = it->substr(5,it->length() - 5);

            if (type_string.compare("DELETION") == 0) {
                type = SV_BreakPoint::DELETION;
            } else if (type_string.compare("DUPLICATION") == 0) {
                type = SV_BreakPoint::DUPLICATION;
            } else if (type_string.compare("INVERSION") == 0) {
                type = SV_BreakPoint::INVERSION;
            } else if (type_string.compare("TRANSLOCATION") == 0) {
                type = SV_BreakPoint::TRANSLOCATION;
            } else {
                cerr << "ERROR IN BEDPE FILE.  TYPE \""<< type_string <<
                     "\" not supported " <<
                     "(DELETION,DUPLICATION,INVERSION,TRANSLOCATION)" <<
                     endl;
                abort();
            }

        }
    }

    if (type == -1) {
        cerr << "ERROR IN BEDPE FILE.  Either no TYPE field. " <<
             endl;
        abort();
    }

    weight = _weight;
}
//}}}

//{{{ ostream& operator << (ostream& out, const SV_Pair& p)
ostream& operator << (ostream& out, const SV_Bedpe& p)
{

    out << p.side_l.chr << "," <<
        p.side_l.start << "," <<
        p.side_l.end << "," <<
        p.side_l.strand <<
        "\t" <<
        p.side_r.chr << "," <<
        p.side_r.start << "," <<
        p.side_r.end << "," <<
        p.side_r.strand;

    return out;
}
//}}}

//{{{ SV_BreakPoint* SV_Pair:: get_bp()
SV_BreakPoint*
SV_Bedpe::
get_bp()
{
    // Make a new break point
    SV_BreakPoint *new_bp = new SV_BreakPoint(this);

    set_bp_interval_start_end(&(new_bp->interval_l),
                              &side_l,
                              &side_r,
                              0,
                              reader->distro_size);
    set_bp_interval_start_end(&(new_bp->interval_r),
                              &side_r,
                              &side_l,
                              0,
                              reader->distro_size);

    //set_bp_interval_probability(&(new_bp->interval_l));
    //set_bp_interval_probability(&(new_bp->interval_r));
    new_bp->interval_r.p = NULL;
    new_bp->interval_l.p = NULL;

    new_bp->type = type;
    new_bp->weight = weight;

    return new_bp;
}
//}}}

//{{{ void SV_Bedpe:: set_bp_interval_start_end(struct breakpoint_interval *i,
void
SV_Bedpe::
set_bp_interval_start_end(struct breakpoint_interval *i,
                          struct interval *target_interval,
                          struct interval *target_pair,
                          int back_distance,
                          int distro_size)
{
    i->i.chr = target_interval->chr;
    i->i.strand = target_interval->strand;
    i->i.start = target_interval->start;
    i->i.end = target_interval->end;
    i->i.start_clip = 0;
}
//}}}

//{{{ log_space* SV_Bedpe:: get_bp_interval_probability(char strand)
log_space*
SV_Bedpe::
get_bp_interval_probability(char strand,
                            int distro_size,
                            double *distro)
{
    CHR_POS size = distro_size;
    log_space *tmp_p = (log_space *) malloc(size * sizeof(log_space));
    CHR_POS j;
    for (j = 0; j < size; ++j) {
        if (strand == '+')
            tmp_p[j] = get_ls(distro[j]);
        else
            tmp_p[(size - 1) - j] = get_ls(distro[j]);
    }

    return tmp_p;
}
//}}}

//{{{ void SV_Pair:: print_evidence()
void
SV_Bedpe::
print_evidence()
{
    print_bedpe(0);
}
//}}}

//{{{ void SV_Bedpe:: print_bedpe(int score)
void
SV_Bedpe::
print_bedpe(int score)
{
    // use the address of the current object as the id
    string sep = "\t";
    cout <<
         side_l.chr << sep <<
         side_l.start << sep <<
         (side_l.end + 1) << sep <<
         side_r.chr << sep <<
         side_r.start << sep<<
         (side_r.end + 1) << sep<<
         this << sep <<
         score << sep <<
         side_l.strand << "\t" <<
         side_r.strand <<
         endl;
}
//}}}

//{{{ void SV_Bedpe:: process_bedpe(const BEDPE *bedpeEntry,
void
SV_Bedpe::
process_bedpe(const BEDPE *bedpeEntry,
              UCSCBins<SV_BreakPoint*> &r_bin,
              int weight,
              int ev_id,
              SV_BedpeReader *reader)
{
    SV_Bedpe *new_bedpe = new SV_Bedpe(bedpeEntry,
                                       weight,
                                       ev_id,
                                       reader);

    SV_BreakPoint *new_bp = new_bedpe->get_bp();
    new_bp->cluster(r_bin);
}
//}}}

//{{{ string SV_Bedpe:: evidence_type()
string
SV_Bedpe::
evidence_type()
{
    return "bedpe";
}
//}}}

//{{{ void SV_Pair:: set_interval_probability()
void
SV_Bedpe::
set_bp_interval_probability(struct breakpoint_interval *i)
{
    CHR_POS size = i->i.end - i->i.start + 1;
    log_space *tmp_p = (log_space *) malloc(size * sizeof(log_space));
    //log_space *src_p;

    double v = 1.0 / ((double )size);

    //int offset = i->i.start_clip;
    CHR_POS j;
    for (j = 0; j < size; ++j) {
        //tmp_p[j] = src_p[j + offset];
        tmp_p[j] = get_ls(v);
    }

    i->p = tmp_p;
}
//}}}
