/*****************************************************************************
 * SV_Pair.cpp
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
#include "api/BamWriter.h"
using namespace BamTools;

#include "SV_PairReader.h"
#include "SV_BreakPoint.h"
#include "SV_Pair.h"
#include "SV_Tools.h"
#include "log_space.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <math.h>

using namespace std;

//{{{ statics
double  SV_Pair:: insert_Z = 0;
//}}}

//{{{ SV_Pair:: SV_Pair(const BamAlignment &bam_a,
// if both reads are on the same chrome, then read_l must map before read_r
// if the reads are on different strands then read_l must be on the lexo
// lesser chrom (using the string.compare() method)
SV_Pair::
SV_Pair(const BamAlignment &bam_a,
        const BamAlignment &bam_b,
        const RefVector &refs,
        int _weight,
        int _ev_id,
        SV_PairReader *_reader)
{
    reader = _reader;

    if ( bam_a.MapQuality < bam_b.MapQuality )
        min_mapping_quality = bam_a.MapQuality;
    else
        min_mapping_quality = bam_b.MapQuality;

    struct interval tmp_a, tmp_b;
    tmp_a.start = bam_a.Position;
    tmp_a.end = bam_a.GetEndPosition(false, false) - 1;
    tmp_a.chr = refs.at(bam_a.RefID).RefName;

    if ( bam_a.IsReverseStrand() == true )
        tmp_a.strand = '-';
    else
        tmp_a.strand = '+';

    tmp_b.start = bam_b.Position;
    tmp_b.end = bam_b.GetEndPosition(false, false) - 1;
    tmp_b.chr = refs.at(bam_b.RefID).RefName;

    if ( bam_b.IsReverseStrand() == true )
        tmp_b.strand = '-';
    else
        tmp_b.strand = '+';

    //if ( tmp_a.chr.compare(tmp_b.chr) > 0 ) {
    if ( bam_a.RefID < bam_b.RefID ) {
        read_l = tmp_a;
        read_r = tmp_b;
        //} else if ( tmp_a.chr.compare(tmp_b.chr) < 0 ) {
    } else if ( bam_a.RefID > bam_b.RefID) {
        read_l = tmp_b;
        read_r = tmp_a;
    } else { // ==
        if (tmp_a.start > tmp_b.start) {
            read_l = tmp_b;
            read_r = tmp_a;
        } else {
            read_l = tmp_a;
            read_r = tmp_b;
        }
    }

    weight = _weight;
    ev_id = _ev_id;
}
//}}}

//{{{ log_space* SV_Pair:: get_interval_probability(char strand)
log_space*
SV_Pair::
get_bp_interval_probability(char strand,
                            int distro_size,
                            double *distro)
{
    int size = distro_size;
    log_space *tmp_p = (log_space *) malloc(size * sizeof(log_space));
    int j;
    for (j = 0; j < size; ++j) {
        if (strand == '+')
            tmp_p[j] = get_ls(distro[j]);
        else
            tmp_p[(size - 1) - j] = get_ls(distro[j]);
    }

    return tmp_p;
}
//}}}

//{{{ void SV_Pair:: set_bp_interval_start_end(struct breakpoint_interval *i,
/* targer_pair is the other interval in the pair, it is used to bound the
 * breakpoint probabilty distribution range.  The breakpoint cannot exist after
 * the start of the targer_pair
 */
void
SV_Pair::
set_bp_interval_start_end(struct breakpoint_interval *i,
                          struct interval *target_interval,
                          struct interval *target_pair,
                          unsigned int back_distance,
                          unsigned int distro_size)
{
    i->i.chr = target_interval->chr;
    i->i.strand = target_interval->strand;
    if ( i->i.strand == '+' ) {

        if (back_distance > target_interval->end) {
            i->i.start = 0;
            i->i.start_clip = target_interval->end;
        }
        else {
            i->i.start = target_interval->end - back_distance;
            i->i.start_clip = 0;
        }

        i->i.end = i->i.start + distro_size - 1;

    } else {
        i->i.end = target_interval->start + back_distance;

        if (distro_size > i->i.end) {
            i->i.start = 0;
            i->i.start_clip = i->i.end;
        }
        else {
            i->i.start = i->i.end - distro_size + 1;
            i->i.start_clip = 0;
        }
    }
}
//}}}

//{{{ SV_BreakPoint* SV_Pair:: get_bp()
SV_BreakPoint*
SV_Pair::
get_bp()
{
    // Make a new break point
    SV_BreakPoint *new_bp = new SV_BreakPoint(this);

    set_bp_interval_start_end(&(new_bp->interval_l),
                              &read_l,
                              &read_r,
                              reader->back_distance,
                              reader->distro_size);
    set_bp_interval_start_end(&(new_bp->interval_r),
                              &read_r,
                              &read_l,
                              reader->back_distance,
                              reader->distro_size);

    // test to see if the two interval intersect the reads
    if ( (new_bp->interval_l.i.chr.compare(read_r.chr) == 0 ) &&
            (new_bp->interval_l.i.end >= read_r.start) ) {
        new_bp->interval_l.i.end_clip =
            new_bp->interval_l.i.end - read_r.start + 1;
        new_bp->interval_l.i.end = read_r.start - 1;
    }

    if ( (new_bp->interval_r.i.chr.compare(read_l.chr) == 0 ) &&
            (new_bp->interval_r.i.start <= read_l.end) ) {
        new_bp->interval_r.i.start_clip =
            read_l.end - new_bp->interval_r.i.start + 1;
        new_bp->interval_r.i.start = read_l.end + 1;
    }

    new_bp->interval_r.p = NULL;
    new_bp->interval_l.p = NULL;

    if (new_bp->interval_l.i.strand == '+') // // + ?
        if (new_bp->interval_r.i.strand == '+') // + +
            new_bp->type = SV_BreakPoint::INVERSION;
        else // + -
            new_bp->type = SV_BreakPoint::DELETION;
    else // - ?
        if (new_bp->interval_r.i.strand == '+') // - +
            new_bp->type = SV_BreakPoint::DUPLICATION;
        else // - -
            new_bp->type = SV_BreakPoint::INVERSION;

    new_bp->weight = weight;

    return new_bp;
}
//}}}

//{{{ bool SV_Pair:: is_aberrant()
bool
SV_Pair::
is_aberrant()
{
    if ( read_l.strand == read_r.strand )
        return true;

    if ( read_l.strand == '-')
        return true;

    if ( (read_r.end - read_l.start) >=
            reader->mean + (reader->discordant_z*reader->stdev) )
        return true;

    if ( (read_r.end - read_l.start) <=
            reader->mean - (reader->discordant_z*reader->stdev) )
        return true;

    return false;
}
//}}}

//{{{ bool SV_Pair:: is_sane()
bool
SV_Pair::
is_sane()
{
    // test if either end is in and excluded region
    vector<UCSCElement<int> > v = exclude_regions.get(read_l.chr,
                                  read_l.start,
                                  read_l.end,
                                  '+',
                                  false);
    if ( v.size() > 0 )
        return false;

    v = exclude_regions.get(read_r.chr,
                            read_r.start,
                            read_r.end,
                            '+',
                            false);

    if ( v.size() > 0 )
        return false;


    if ( min_mapping_quality < reader->min_mapping_threshold )
        return false;

    int read_len_a = read_l.end - read_l.start;
    int read_len_b = read_r.end - read_r.start;

    // |-----------------|
    //              |-----------------|
    //              ^min ^max
    //				|-----| overlap
    // |------------| non_overlap
    //
    // |-----------------|
    //                            |-----------------|
    //                   ^min     ^max
    //                   |--------| "overalp" (less than 0)
    // |--------| "non_overlap"
    //

#ifdef TRACE
    cerr << "READ LEN\t" <<
            read_len_a << "\t" <<
            read_len_b << endl;
#endif 

    if ( (read_len_a <= 0) || (read_len_b <= 0) )
        return false;


    int overlap = min(read_l.end, read_r.end) - max(read_l.start, read_r.start);
    //int non_overlap = min(read_len_a, read_len_b) - overlap;
    // how much does not overlap // overlap is at most read_len
    //
#ifdef TRACE
    cerr << "OVERLAP\t" <<
            read_l.chr << " " << 
            read_l.start  << " " << 
            read_l.end  << "\t" << 
            read_r.chr << " " << 
            read_r.start  << " " << 
            read_r.end  << "\t" <<
            overlap << endl;
#endif 

    if (overlap > 0)
        return false;
    else
        return true;
    /*
    int non_overlap = min(read_len_a, read_len_b) - overlap;

    if ( (overlap > 0) && (abs(non_overlap) < reader->min_non_overlap) )
        return false;
    else
        return true;
    */
}
//}}}

//{{{ ostream& operator << (ostream& out, const SV_Pair& p)
ostream& operator << (ostream& out, const SV_Pair& p)
{

    out << p.read_l.chr << "," <<
        p.read_l.start << "," <<
        p.read_l.end << "," <<
        p.read_l.strand <<
        "\t" <<
        p.read_r.chr << "," <<
        p.read_r.start << "," <<
        p.read_r.end << "," <<
        p.read_r.strand;

    return out;
}
//}}}

//{{{ void SV_Pair:: print_evidence()
void
SV_Pair::
print_evidence()
{
    print_bedpe(0);
}
//}}}

//{{{ void BreakPoint:: print_bedpe()
void
SV_Pair::
print_bedpe(int score)
{
    // use the address of the current object as the id
    string sep = "\t";
    cout <<
         read_l.chr << sep <<
         read_l.start << sep <<
         read_l.end << sep <<
         read_r.chr << sep <<
         read_r.start << sep<<
         read_r.end << sep<<
         this << sep <<
         score << sep <<
         read_l.strand << sep <<
         read_r.strand << sep <<
	 "id:" << ev_id << sep <<
         "weight:" << weight <<
         endl;
}
//}}}

//{{{ set_distro_from_histo (int back_distance,
// We will assume that the distrobution to upsstream will follow the histogram
// and downstream will be follow an exponetial decay distribution based on the
// back_distance
int
SV_Pair::
set_distro_from_histo (int back_distance,
                       int histo_start,
                       int histo_end,
                       double *histo,
                       double **distro)
{
    // It is possible for the histo to contain 0 values at the ends,
    // to prevent further issues, we will trim them here
    /*
    for (int i = 0; i < histo_end - histo_start + 1; ++i) {
        if (histo[i] == 0) {
            histo_end = i-1;
            break;
        }
    }
    */

    double lambda = log(0.0001)/(-1 * back_distance);
    // the bp distribution begins SV_Pair::back_distance base pairs back
    // from the end of the read (or begining for the negative strand), then
    // extends for SV_Pair::histo_end base pairs, so the total size is
    // SV_Pair::back_distance + SV_Pair::histo_end
    int distro_size = back_distance + histo_end;

    *distro = (double *) malloc(distro_size * sizeof(double));

    for (int i = 0; i < back_distance; ++i)
        (*distro)[i] = exp(-1*lambda*(back_distance - i));

    for (int i = back_distance; i < histo_start; ++i)
        (*distro)[i] = 1;

    double last = 0;
    for (int i = histo_end - 1; i >= histo_start; --i) {
        (*distro)[i + back_distance] = histo[i - histo_start] + last;
        last = (*distro)[i + back_distance];
    }


    int zero_i = -1;
    for (int i = 0; i < distro_size; ++i) {
        if ((*distro)[i] == 0) {
            cerr << i << "\t" << (*distro)[i] << endl;
            zero_i = i;
            break;
        }
    }

    if (zero_i != -1) {
        distro_size = zero_i;
        double *new_distro = (double *) malloc(distro_size * sizeof(double));
        for (int i = 0; i < distro_size; ++i) {
            new_distro[i] = (*distro)[i];
        }
        free(*distro);
        *distro = new_distro;
    }


    return distro_size;
}
//}}}

//{{{ void process_pair(const BamAlignment &curr,
void
SV_Pair::
process_pair(const BamAlignment &curr,
             const RefVector refs,
             map<string, BamAlignment> &mapped_pairs,
             UCSCBins<SV_BreakPoint*> &r_bin,
             int weight,
             int ev_id,
             SV_PairReader *reader)
{
    if (mapped_pairs.find(curr.Name) == mapped_pairs.end())
        mapped_pairs[curr.Name] = curr;
    else {
        SV_Pair *new_pair = new SV_Pair(mapped_pairs[curr.Name],
                                        curr,
                                        refs,
                                        weight,
                                        ev_id,
                                        reader);
        //cerr << count_clipped(curr.CigarData) << "\t" <<
                //count_clipped(mapped_pairs[curr.Name].CigarData) << endl;
                
        if ( new_pair->is_sane() &&  
             new_pair->is_aberrant() &&
             (count_clipped(curr.CigarData) > 0) &&
             (count_clipped(mapped_pairs[curr.Name].CigarData) > 0) ) {
            SV_BreakPoint *new_bp = new_pair->get_bp();
#ifdef TRACE
            cerr << "READ\t" << 
                    refs.at(mapped_pairs[curr.Name].RefID).RefName << "," <<
                    mapped_pairs[curr.Name].Position << "," <<
                    (mapped_pairs[curr.Name].GetEndPosition(false, false) - 1)
                        << "\t" <<
                    refs.at(curr.RefID).RefName << "," <<
                    curr.Position << "," <<
                    (curr.GetEndPosition(false, false) - 1)
                        <<
                    endl;

            cerr << "\tPE\t" << *new_bp << endl;
#endif
            new_bp->cluster(r_bin);
        } else {
            delete(new_pair);
        }
        mapped_pairs.erase(curr.Name);
    }
}
//}}}

//{{{ void process_intra_chrom_pair(const BamAlignment &curr,
void
SV_Pair::
process_intra_chrom_pair(const BamAlignment &curr,
                         const RefVector refs,
                         BamWriter &inter_chrom_reads,
                         map<string, BamAlignment> &mapped_pairs,
                         UCSCBins<SV_BreakPoint*> &r_bin,
                         int weight,
                         int ev_id,
                         SV_PairReader *reader)
{
    if (curr.RefID == curr.MateRefID) {

        process_pair(curr,
                     refs,
                     mapped_pairs,
                     r_bin,
                     weight,
                     ev_id,
                     reader);

    } else if (curr.IsMapped() &&
               curr.IsMateMapped() &&
               (curr.RefID >= 0) &&
               (curr.MateRefID >= 0) ) {
        BamAlignment al = curr;
        string x = reader->get_source_file_name();
        al.AddTag("LS","Z",x);
        inter_chrom_reads.SaveAlignment(al);
    }
}
//}}}

//{{{string SV_Pair:: evidence_type()
string
SV_Pair::
evidence_type()
{
    return "Pair";
}
//}}}

#if 0
//{{{ void SV_Pair:: update_matrix(boost::numeric::ublas::matrix<log_space> *m)
void
SV_Pair::
update_matrix(boost::numeric::ublas::matrix<log_space> *m)
{
    for (unsigned int a = 0; a < m->size1 (); ++ a) {
        for (unsigned int b = 0; b < m->size2 (); ++ b) {
            unsigned int test_break_a, test_break_b, a_d, b_d;

            // Get the distance from left tag to left break and break to
            // right tag, then test the likelyhood of that distance given
            // the known distribution (mean, stdev)
            test_break_a = read_l.start + a;
            test_break_b = read_r.start + b;

            if (read_l.strand == '+')
                a_d = test_break_a - read_l.start;
            else
                a_d = read_l.end - test_break_a;


            if (read_r.strand == '+')
                b_d = test_break_b - read_r.start;
            else
                b_d = read_r.end - test_break_b;

            double dist = (a_d + b_d);
            double z_dist = ((a_d + b_d) - insert_mean);

            double p = gsl_ran_gaussian_pdf(z_dist, insert_stdev);
            log_space lp = get_ls(p);

            (*m)(a, b) = ls_multiply((*m)(a, b), lp);
        }
    }
}
//}}}
//{{{void set_sv_pair_distro()
// We will assume that the distrobution to upsstream will follow the histogram
// and downstream will be follow an exponetial decay distribution based on the
// back_distance
void
SV_Pair::
set_distro_from_histo ()
{
    double lambda = log(0.0001)/(-1 * SV_Pair::back_distance);
    // the bp distribution begins SV_Pair::back_distance base pairs back
    // from the end of the read (or begining for the negative strand), then
    // extends for SV_Pair::histo_end base pairs, so the total size is
    // SV_Pair::back_distance + SV_Pair::histo_end
    SV_Pair::distro_size = SV_Pair::back_distance + SV_Pair::histo_end;
    //SV_Pair::distro = (double *)
    SV_Pair::distro = (double *)
                      malloc(SV_Pair::distro_size * sizeof(double));

    for (int i = 0; i < SV_Pair::back_distance; ++i)
        SV_Pair::distro[i] = exp(-1*lambda*(SV_Pair::back_distance - i));

    for (int i = SV_Pair::back_distance; i < SV_Pair::histo_start; ++i)
        SV_Pair::distro[i] = 1;

    double last = 0;
    for (int i = SV_Pair::histo_end - 1; i >= SV_Pair::histo_start; --i) {
        SV_Pair::distro[i + SV_Pair::back_distance] =
            SV_Pair::histo[i - SV_Pair::histo_start] + last;
        last = SV_Pair::distro[i + SV_Pair::back_distance];
    }
}
//}}}
#endif


