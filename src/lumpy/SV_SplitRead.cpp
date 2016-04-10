/*****************************************************************************
 * SV_SplitRead.cpp
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

#include <exception>
#include "SV_SplitReadReader.h"
#include "SV_BreakPoint.h"
#include "SV_SplitRead.h"
#include "SV_Tools.h"
#include "log_space.h"
#include "bedFilePE.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <math.h>

using namespace std;

//{{{int SV_SplitRead:: count_clipped(vector< CigarOp > cigar_data)
uint32_t
SV_SplitRead::
count_clipped(vector< CigarOp > cigar_data)
{
    uint32_t clipped_count = 0;
    vector< CigarOp >::iterator i;
    for (i = cigar_data.begin(); i != cigar_data.end(); ++i) {
        if ( (i->Type == 'H') || (i->Type == 'S') ) {
            uint32_t l = i->Length;
            clipped_count += l;
        }
    }

    return clipped_count;
}
//}}}

//{{{ void SV_SplitRead:: update_cigar_query(CigarOp op,
void
SV_SplitRead::
update_cigar_query(CigarOp op,
                   int *op_position,
                   struct cigar_query *q)
{
    //cerr << op.Type << op.Length << endl;
    if ( ( *op_position == 0 ) &&
            ( ( op.Type == 'H' ) || ( op.Type == 'S' ) ) ) {
        q->qs_pos = q->qs_pos + op.Length;
        q->qe_pos = q->qe_pos + op.Length;
        q->q_len = q->q_len + op.Length;
    } else if ( ( op.Type == 'H' ) || ( op.Type == 'S' ) ) {
        q->q_len = q->q_len + op.Length;
    } else if ( ( op.Type == 'M' ) || ( op.Type == 'I' ) ) {
        q->qe_pos = q->qe_pos + op.Length;
        q->q_len = q->q_len + op.Length;
    }

    *op_position = *op_position + 1;
}
//}}}

//{{{struct cigar_query SV_SplitRead:: calc_query_pos_from_cigar(vector<
struct cigar_query
        SV_SplitRead::
calc_query_pos_from_cigar(vector< CigarOp > cigar_data,
                          bool is_reverse_strand)
{
    struct cigar_query new_query;
    new_query.qs_pos = 0;
    new_query.qe_pos = 0;
    new_query.q_len = 0;
    int op_position = 0;


    if (is_reverse_strand == true) {
        vector< CigarOp >::reverse_iterator i;
        for (i = cigar_data.rbegin(); i != cigar_data.rend(); ++i)
            update_cigar_query(*i, &op_position, &new_query);
    } else {
        vector< CigarOp >::iterator i;
        for (i = cigar_data.begin(); i != cigar_data.end(); ++i)
            update_cigar_query(*i, &op_position, &new_query);
    }

    return new_query;
}
//}}}

//{{{ SV_SplitRead:: SV_SplitRead(vector< BamAlignment > &block,
SV_SplitRead::
SV_SplitRead(vector< BamAlignment > &block,
             const RefVector &refs,
             int _weight,
             int _ev_id,
             SV_SplitReadReader *_reader)
{
    reader = _reader;
    ev_id = _ev_id;

    if ( block.at(0).MapQuality < block.at(1).MapQuality )
        min_mapping_quality = block.at(0).MapQuality;
    else
        min_mapping_quality = block.at(1).MapQuality;



    struct cigar_query query_a =
        calc_query_pos_from_cigar(block.at(0).CigarData,
                                  block.at(0).IsReverseStrand() );
    struct cigar_query query_b =
        calc_query_pos_from_cigar(block.at(1).CigarData,
                                  block.at(1).IsReverseStrand() );

    struct interval tmp_a, tmp_b;

    tmp_a.strand = '+';
    if (block.at(0).IsReverseStrand())
        tmp_a.strand = '-';
    tmp_a.chr = refs.at(block.at(0).RefID).RefName;
    tmp_a.start = block.at(0).Position;
    tmp_a.end = block.at(0).GetEndPosition();

    tmp_b.strand = '+';
    if (block.at(1).IsReverseStrand())
        tmp_b.strand = '-';
    tmp_b.chr = refs.at(block.at(1).RefID).RefName;
    tmp_b.start = block.at(1).Position;
    tmp_b.end = block.at(1).GetEndPosition();


    if ( (block.at(0).RefID > block.at(1).RefID) ||
            ( (block.at(0).RefID == block.at(1).RefID) &&
              (tmp_a.start > tmp_b.start ) ) ) {
        side_r = tmp_a;
        side_l = tmp_b;
        query_r = query_a;
        query_l = query_b;
    } else {
        side_l = tmp_a;
        side_r = tmp_b;
        query_l = query_a;
        query_r = query_b;
    }

    if (side_l.strand != side_r.strand)
        type = SV_BreakPoint::INVERSION;
    else if ( (	( side_l.strand == '+' ) &&
                ( side_r.strand == '+' ) &&
                ( query_l.qs_pos < query_r.qs_pos ) ) ||
              (	( side_l.strand == '-' ) &&
                  ( side_r.strand == '-' ) &&
                  ( query_l.qs_pos > query_r.qs_pos) ) )
        type = SV_BreakPoint::DELETION;
    else if ( ( ( side_l.strand == '+' ) &&
                ( side_r.strand == '+' ) &&
                ( query_l.qs_pos > query_r.qs_pos ) ) ||
              ( ( side_l.strand == '-' ) &&
                ( side_r.strand == '-' ) &&
                ( query_l.qs_pos < query_r.qs_pos) ) )
        type = SV_BreakPoint::DUPLICATION;
    else {
        cerr << "ERROR IN BAM FILE.  " <<
             "TYPE not detected (DELETION,DUPLICATION,INVERSION)" <<
             endl;
        cerr << "\t" << query_l.qs_pos << "," << side_l.strand << "\t" <<
             query_r.qs_pos << "," << side_r.strand << "\t" <<
             tmp_a.chr << "," << tmp_a.start << "," << tmp_a.end << "\t" <<
             tmp_b.chr << "," << tmp_b.start << "," << tmp_b.end << "\t" <<
             endl;

        throw(1);
    }

    weight = _weight;
}
//}}}

//{{{ SV_SplitRead:: SV_SplitRead(vector< BamAlignment > &block,
SV_SplitRead::
SV_SplitRead(const BamAlignment &bam_a,
             const BamAlignment &bam_b,
             const RefVector &refs,
             int _weight,
             int _ev_id,
             SV_SplitReadReader *_reader)
{
    reader = _reader;
    ev_id = _ev_id;

    // Add readId information to track read evidence. Note that both bam_a and
    // bam_b refer to the same read, thus have the same readId
    read_id = bam_a.Name;

    if ( bam_a.MapQuality < bam_b.MapQuality )
        min_mapping_quality = bam_a.MapQuality;
    else
        min_mapping_quality = bam_b.MapQuality;

    struct cigar_query query_a =
        calc_query_pos_from_cigar(bam_a.CigarData,
                                  bam_a.IsReverseStrand() );
    struct cigar_query query_b =
        calc_query_pos_from_cigar(bam_b.CigarData,
                                  bam_b.IsReverseStrand() );

    struct interval tmp_a, tmp_b;

    tmp_a.strand = '+';
    if (bam_a.IsReverseStrand())
        tmp_a.strand = '-';
    tmp_a.chr = refs.at(bam_a.RefID).RefName;
    tmp_a.start = bam_a.Position;
    tmp_a.end = bam_a.GetEndPosition();

    tmp_b.strand = '+';
    if (bam_b.IsReverseStrand())
        tmp_b.strand = '-';
    tmp_b.chr = refs.at(bam_b.RefID).RefName;
    tmp_b.start = bam_b.Position;
    tmp_b.end = bam_b.GetEndPosition();


    //if ( ( tmp_a.chr.compare(tmp_b.chr) > 0 ) ||
    //( ( tmp_a.chr.compare(tmp_b.chr) == 0 ) &&
    //( tmp_a.start > tmp_b.start ) ) ) {

    if ( (bam_a.RefID > bam_b.RefID) ||
            ( (bam_a.RefID == bam_b.RefID) &&
              (tmp_a.start > tmp_b.start ) ) ) {
        side_r = tmp_a;
        side_l = tmp_b;
        query_r = query_a;
        query_l = query_b;
    } else {
        side_l = tmp_a;
        side_r = tmp_b;
        query_l = query_a;
        query_r = query_b;
    }

    if (side_l.strand != side_r.strand)
        type = SV_BreakPoint::INVERSION;
    else if ( (	( side_l.strand == '+' ) &&
                ( side_r.strand == '+' ) &&
                ( query_l.qs_pos < query_r.qs_pos ) ) ||
              (	( side_l.strand == '-' ) &&
                  ( side_r.strand == '-' ) &&
                  ( query_l.qs_pos > query_r.qs_pos) ) )
        type = SV_BreakPoint::DELETION;
    else if ( ( ( side_l.strand == '+' ) &&
                ( side_r.strand == '+' ) &&
                ( query_l.qs_pos > query_r.qs_pos ) ) ||
              ( ( side_l.strand == '-' ) &&
                ( side_r.strand == '-' ) &&
                ( query_l.qs_pos < query_r.qs_pos) ) )
        type = SV_BreakPoint::DUPLICATION;
    else {
        cerr << "ERROR IN BAM FILE.  " <<
             "TYPE not detected (DELETION,DUPLICATION,INVERSION)" <<
             endl;
        cerr << "\t" << query_l.qs_pos << "," << side_l.strand << "\t" <<
             query_r.qs_pos << "," << side_r.strand << "\t" <<
             tmp_a.chr << "," << tmp_a.start << "," << tmp_a.end << "\t" <<
             tmp_b.chr << "," << tmp_b.start << "," << tmp_b.end << "\t" <<
             endl;

        throw(1);
    }

    weight = _weight;
}
//}}}

//{{{ SV_BreakPoint* SV_SplitRead:: get_bp()
SV_BreakPoint*
SV_SplitRead::
get_bp()
{
    // Make a new break point
    SV_BreakPoint *new_bp = new SV_BreakPoint(this);

    new_bp->type = type;

    new_bp->interval_l.i.chr = side_l.chr;
    new_bp->interval_l.i.strand = '?';

    new_bp->interval_r.i.chr = side_r.chr;
    new_bp->interval_r.i.strand = '?';

    if (type == SV_BreakPoint::INVERSION) {
        if ( ( ( query_l.qs_pos < query_r.qs_pos ) &&
                ( side_l.strand == '-' ) &&
                ( side_r.strand == '+' ) ) ||
                ( ( query_l.qs_pos > query_r.qs_pos ) &&
                  ( side_l.strand == '+' ) &&
                  ( side_r.strand == '-' ) ) ) {

            new_bp->interval_l.i.strand = '-';

            if ( side_l.start > reader->back_distance) {
                new_bp->interval_l.i.start =
                    side_l.start - reader->back_distance;
                new_bp->interval_l.i.start_clip = 0;
            } else {
                new_bp->interval_l.i.start = 0;
                new_bp->interval_l.i.start_clip = side_l.start;
            }

            new_bp->interval_l.i.end = side_l.start + reader->back_distance;


            new_bp->interval_r.i.strand = '-';

            if (side_r.start > reader->back_distance) {
                new_bp->interval_r.i.start =
                    side_r.start - reader->back_distance;
                new_bp->interval_r.i.start_clip = 0;
            } else {
                new_bp->interval_r.i.start = 0;
                new_bp->interval_r.i.start_clip = side_r.start;
            }

            new_bp->interval_r.i.end = side_r.start + reader->back_distance;

        } else if ( ( ( query_l.qs_pos < query_r.qs_pos ) &&
                      ( side_l.strand == '+' ) &&
                      ( side_r.strand == '-' ) ) ||
                    ( ( query_l.qs_pos > query_r.qs_pos) &&
                      ( side_l.strand == '-' ) &&
                      ( side_r.strand == '+' ) ) ) {

            new_bp->interval_l.i.strand = '+';

            if (side_l.end > reader->back_distance) {
                new_bp->interval_l.i.start =
                    side_l.end  - reader->back_distance;
                new_bp->interval_l.i.start_clip = 0;
            } else {
                new_bp->interval_l.i.start = 0;
                new_bp->interval_l.i.start_clip = side_l.end;
            }

            new_bp->interval_l.i.end = side_l.end + reader->back_distance;

            new_bp->interval_r.i.strand = '+';
            if (side_r.end > reader->back_distance) {
                new_bp->interval_r.i.start =
                    side_r.end - reader->back_distance;
                new_bp->interval_r.i.start_clip = 0;
            } else {
                new_bp->interval_r.i.start = 0;
                new_bp->interval_r.i.start_clip = side_r.end;
            }

            new_bp->interval_r.i.end = side_r.end + reader->back_distance;
        } else {
            cerr << "Cannot determine BP type:" <<
                    "side_l.start:" << side_l.start << "\t" <<
                    "side_l.end:" << side_l.end << "\t" <<
                    "side_r.start:" << side_r.start << "\t" <<
                    "side_r.end:" << side_r.end << "\t" <<
                    "query_l.qs_pos:" << query_l.qs_pos << "\t" <<
                    "query_r.qs_pos:" << query_r.qs_pos<< "\t" <<
                    "side_l.strand:" << side_l.strand<< "\t" <<
                    "side_r.strand:" << side_r.strand<< endl;
            return NULL;
            //abort();
        }
    } else if (type == SV_BreakPoint::DELETION) {

        new_bp->interval_l.i.strand = '+';

        if (side_l.end > reader->back_distance) {
            new_bp->interval_l.i.start = side_l.end - reader->back_distance;
            new_bp->interval_l.i.start_clip = 0;
        } else {
            new_bp->interval_l.i.start = 0;
            new_bp->interval_l.i.start_clip = side_l.end;
        }

        new_bp->interval_l.i.end = side_l.end + reader->back_distance;

        new_bp->interval_r.i.strand = '-';

        if (side_r.start > reader->back_distance) {
            new_bp->interval_r.i.start = side_r.start - reader->back_distance;
            new_bp->interval_r.i.start_clip = 0;
        } else {
            new_bp->interval_r.i.start = 0;
            new_bp->interval_r.i.start_clip = side_r.start;
        }

        new_bp->interval_r.i.end = side_r.start + reader->back_distance;
    } else if (type == SV_BreakPoint::DUPLICATION) {

        new_bp->interval_l.i.strand = '-';

        if (side_l.start > reader->back_distance) {
            new_bp->interval_l.i.start = side_l.start - reader->back_distance;
            new_bp->interval_l.i.start_clip = 0;
        } else {
            new_bp->interval_l.i.start = 0;
            new_bp->interval_l.i.start_clip = side_l.start;
        }

        new_bp->interval_l.i.end = side_l.start + reader->back_distance;

        new_bp->interval_r.i.strand = '+';

        if (side_r.end > reader->back_distance) {
            new_bp->interval_r.i.start = side_r.end - reader->back_distance;
            new_bp->interval_r.i.start_clip = 0;
        } else {
            new_bp->interval_r.i.start = 0;
            new_bp->interval_r.i.start_clip = side_r.end;
        }

        new_bp->interval_r.i.end = side_r.end + reader->back_distance;
    }  else {
        cerr << "Cannot determine BP type:" <<
                "side_l.start:" << side_l.start << "\t" <<
                "side_l.end:" << side_l.end << "\t" <<
                "side_r.start:" << side_r.start << "\t" <<
                "side_r.end:" << side_r.end << "\t" <<
                "query_l.qs_pos:" << query_l.qs_pos << "\t" <<
                "query_r.qs_pos:" << query_r.qs_pos<< "\t" <<
                "side_l.strand:" << side_l.strand<< "\t" <<
                "side_r.strand:" << side_r.strand<< endl;
        return NULL;
        //abort();
    }

    if ( (new_bp->interval_l.i.chr.compare(side_r.chr) == 0 ) &&
            (new_bp->interval_l.i.end >= side_r.start) ) {
        new_bp->interval_l.i.end_clip =
            new_bp->interval_l.i.end - side_r.start + 1;
        new_bp->interval_l.i.end = side_r.start - 1;
    }

    if ( (new_bp->interval_r.i.chr.compare(side_l.chr) == 0) &&
            (new_bp->interval_r.i.start <= side_l.end)) {
        new_bp->interval_r.i.start_clip =
            side_l.end - new_bp->interval_r.i.start + 1;
        new_bp->interval_r.i.start = side_l.end + 1;
    }

    /*
    cerr <<
    		side_l.chr << "," <<
    		side_l.start << "," <<
    		side_l.end << "\t" <<
    		side_r.chr << "," <<
    		side_r.start << "," <<
    		side_r.end << "\t" <<

    		new_bp->interval_l.i.chr << "," <<
    		new_bp->interval_l.i.start << "," <<
    		new_bp->interval_l.i.end << "\t" <<
    		new_bp->interval_r.i.chr << "," <<
    		new_bp->interval_r.i.start << "," <<
    		new_bp->interval_r.i.end << "\t" <<
    		(new_bp->interval_l.i.start < new_bp->interval_l.i.end )<< "\t" <<
    		(new_bp->interval_r.i.start < new_bp->interval_r.i.end )<< endl;
    */



    //set_bp_interval_probability(&(new_bp->interval_l));
    //set_bp_interval_probability(&(new_bp->interval_r));
    new_bp->interval_r.p = NULL;
    new_bp->interval_l.p = NULL;

    new_bp->weight = weight;

    return new_bp;
}
//}}}

//{{{ log_space* SV_SplitRead:: get_bp_interval_probability(char strand)
log_space*
SV_SplitRead::
get_bp_interval_probability(char strand,
                            unsigned int back_distance)
{
    double lambda = log(0.0001)/(-1 * (int)back_distance);

    unsigned int distro_size = 2*back_distance + 1;
    log_space *tmp_p = (log_space *) malloc( distro_size * sizeof(log_space));
    unsigned int j;
    for (j = 0; j < back_distance; ++j)
        tmp_p[j] = get_ls( exp(-1*lambda*(back_distance - j)));

    for (; j < distro_size; ++j)
        tmp_p[j] = get_ls(
                       exp(-1*lambda*(back_distance - (distro_size - j - 1))));


    return tmp_p;
}
//}}}

//{{{ bool SV_SplitRead:: is_sane()
bool
SV_SplitRead::
is_sane()
{
    // test if either end is in and excluded region
    vector<UCSCElement<int> > v = exclude_regions.get(side_l.chr,
                                  side_l.start,
                                  side_l.end,
                                  '+',
                                  false);
    if ( v.size() > 0 )
        return false;

    v = exclude_regions.get(side_r.chr,
                            side_r.start,
                            side_r.end,
                            '+',
                            false);

    if ( v.size() > 0 )
        return false;

    if ( min_mapping_quality < reader->min_mapping_threshold )
        return false;

    if (side_l.strand != side_r.strand)
        ;
    else if ( (	( side_l.strand == '+' ) &&
                ( side_r.strand == '+' ) &&
                ( query_l.qs_pos < query_r.qs_pos ) ) ||
              (	( side_l.strand == '-' ) &&
                  ( side_r.strand == '-' ) &&
                  ( query_l.qs_pos > query_r.qs_pos) ) )
        ;
    else if ( ( ( side_l.strand == '+' ) &&
                ( side_r.strand == '+' ) &&
                ( query_l.qs_pos > query_r.qs_pos ) ) ||
              ( ( side_l.strand == '-' ) &&
                ( side_r.strand == '-' ) &&
                ( query_l.qs_pos < query_r.qs_pos) ) )
        ;
    else {
        return false;
    }

    return (side_l.chr != side_r.chr) || (side_l.end < side_r.start);
}
//}}}

//{{{ ostream& operator << (ostream& out, const SV_Pair& p)
ostream& operator << (ostream& out, const SV_SplitRead& p)
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

//{{{ void SV_SplitRead:: print_evidence()
void
SV_SplitRead::
print_evidence()
{
    cout << read_id << "\t";
    print_bedpe(0);
}
//}}}

//{{{ void BreakPoint:: print_bedpe()
void
SV_SplitRead::
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
         side_r.strand << "\t" <<
         "id:" << ev_id << sep <<
         "weight:" << weight <<
         endl;
}
//}}}

//{{{ void process_split(const BamAlignment &curr,
void
SV_SplitRead::
process_split(const BamAlignment &curr,
              const RefVector refs,
              map<string, BamAlignment> &mapped_splits,
              UCSCBins<SV_BreakPoint*> &r_bin,
              int weight,
              int ev_id,
              SV_SplitReadReader *_reader)
{

    if (mapped_splits.find(curr.Name) == mapped_splits.end()) {
        uint32_t clipped = count_clipped(curr.CigarData);
        if (clipped >= _reader->min_clip)
            mapped_splits[curr.Name] = curr;
    } else {
        try {
            SV_SplitRead *new_split_read =
                new SV_SplitRead(mapped_splits[curr.Name],
                                 curr,
                                 refs,
                                 weight,
                                 ev_id,
                                 _reader);

            SV_BreakPoint *new_bp = NULL;
            if (new_split_read->is_sane() &&
                (count_clipped(curr.CigarData) > 0) &&
                (count_clipped(mapped_splits[curr.Name].CigarData) > 0) ) {
                new_bp = new_split_read->get_bp();

                if (new_bp != NULL) {
                    new_bp->cluster(r_bin);
                } else {
                    cerr << "Alignment name:" << curr.Name << endl;
                    free(new_split_read);
                }
            } else {
                //cerr << "not sane" << endl;
                delete(new_split_read);
            }

        } catch (int) {
            cerr << "Error creating split read: " << endl;
        }

        mapped_splits.erase(curr.Name);
    }
}
//}}}

//{{{ void process_intra_chrom_split(const BamAlignment &curr,
void
SV_SplitRead::
process_intra_chrom_split(const BamAlignment &curr,
                          const RefVector refs,
                          BamWriter &inter_chrom_reads,
                          map<string, BamAlignment> &mapped_splits,
                          UCSCBins<SV_BreakPoint*> &r_bin,
                          int weight,
                          int ev_id,
                          SV_SplitReadReader *_reader)
{

    if (mapped_splits.find(curr.Name) == mapped_splits.end()) {
        uint32_t clipped = count_clipped(curr.CigarData);
        if ( curr.HasTag("YP") == true) {
            uint32_t t;
            curr.GetTag("YP", t);
            if (t == 2)
                mapped_splits[curr.Name] = curr;
        }
        else if (clipped >= _reader->min_clip)
            mapped_splits[curr.Name] = curr;
    } else {
        if ( mapped_splits[curr.Name].RefID == curr.RefID ) {
            try {
                SV_SplitRead *new_split_read =
                    new SV_SplitRead(mapped_splits[curr.Name],
                                     curr,
                                     refs,
                                     weight,
                                     ev_id,
                                     _reader);

                SV_BreakPoint *new_bp = NULL;
                if (new_split_read->is_sane()) {
                    new_bp = new_split_read->get_bp();

                    if (new_bp != NULL) {
                        new_bp->cluster(r_bin);
                    } else {
                        cerr << "Alignment name:" << curr.Name << endl;
                        delete(new_split_read);
                    }
                } else {
                    delete(new_split_read);
                }
            } catch (int) {
                cerr << "Error creating split read: " << endl;
            }

        } else {
            BamAlignment al1 = curr;
            BamAlignment al2 = mapped_splits[curr.Name];

            al1.MateRefID = al2.RefID;
            al2.MateRefID = al1.RefID;

            al1.MatePosition = al2.Position;
            al2.MatePosition = al1.Position;

            string x = _reader->get_source_file_name();

            al1.AddTag("LS","Z",x);
            al2.AddTag("LS","Z",x);

            inter_chrom_reads.SaveAlignment(al1);
            inter_chrom_reads.SaveAlignment(al2);
        }
        mapped_splits.erase(curr.Name);
    }
}
//}}}

//{{{string SV_SplitRead:: evidence_type()
string
SV_SplitRead::
evidence_type()
{
    return "split read";
}
//}}}
