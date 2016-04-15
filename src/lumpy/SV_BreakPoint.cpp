/*****************************************************************************
 * SV_BreakPoint.cpp
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


#include "SV_BreakPoint.h"
#include "SV_Evidence.h"
#include "SV_Tools.h"
#include <sstream>
#include <algorithm>
#include <utility>
#include <vector>
#include <climits>
#include <sstream>
#include <map>


double SV_BreakPoint::p_trim_threshold = 1;;
double SV_BreakPoint::p_merge_threshold = 0;
const int SV_BreakPoint::DELETION;
const int SV_BreakPoint::DUPLICATION;
const int SV_BreakPoint::INVERSION;
const int SV_BreakPoint::TRANSLOCATION;


//{{{ SV_BreakPoint:: SV_BreakPoint(SV_Evidence *e)
SV_BreakPoint::
SV_BreakPoint(SV_Evidence *e)
{
    evidence.push_back(e);

    ev_ids[e->ev_id] = 1;
}
//}}}

//{{{ SV_BreakPoint:: SV_BreakPoint(SV_BreakPoint *a, SV_BreakPoint *b)
SV_BreakPoint::
SV_BreakPoint(SV_BreakPoint *a, SV_BreakPoint *b)
{
    a->init_interval_probabilities();
    b->init_interval_probabilities();

    // a may not overlap b, so we need to check
    // after this ar_overlap_intr will point to the interval in a that
    // overlaps the right interval in b, and al_overlap_intr will point to the
    // interval in a that overlaps the left interval in b
    struct breakpoint_interval *al_overlap_intr = NULL,
                                        *ar_overlap_intr = NULL;

    if ( does_intersect(&(a->interval_l), &(b->interval_l), false) &&
            does_intersect(&(a->interval_r), &(b->interval_r), false) ) {
        al_overlap_intr = &(a->interval_l);
        ar_overlap_intr = &(a->interval_r);
    } else if ( does_intersect(&(a->interval_r), &(b->interval_l), false) &&
                does_intersect(&(a->interval_l), &(b->interval_r), false) ) {
        ar_overlap_intr = &(a->interval_l);
        al_overlap_intr = &(a->interval_r);
    } else {
        throw 0;
    }

    if ( (al_overlap_intr != NULL) && (ar_overlap_intr != NULL) ) {
        CHR_POS r_merged_start,
                r_merged_end,
                l_merged_start,
                l_merged_end;

        log_space *r_merged_prob, *l_merged_prob;

        bool l_merges = test_interval_merge(&(b->interval_l),
                                            al_overlap_intr,
                                            &l_merged_start,
                                            &l_merged_end,
                                            &l_merged_prob);


        bool r_merges = test_interval_merge(&(b->interval_r),
                                            ar_overlap_intr,
                                            &r_merged_start,
                                            &r_merged_end,
                                            &r_merged_prob);

        weight = 0;

        if (r_merges && l_merges) {

            interval_l.i.start = l_merged_start;
            interval_l.i.end = l_merged_end;
            interval_l.p = l_merged_prob;

            interval_r.i.start = r_merged_start;
            interval_r.i.end = r_merged_end;
            interval_r.p = r_merged_prob;

            vector<SV_Evidence*>::iterator ev_it;

            for (ev_it = a->evidence.begin(); ev_it < a->evidence.end(); ++ev_it)
                evidence.push_back(*ev_it);

            for (ev_it = b->evidence.begin(); ev_it < b->evidence.end(); ++ev_it)
                evidence.push_back(*ev_it);

            map<int, int>::iterator ev_id_it;

            for (ev_id_it = a->ev_ids.begin(); ev_id_it != a->ev_ids.end(); ++ev_id_it) {
                if ( ev_ids.find( ev_id_it->first ) == ev_ids.end() )
                    ev_ids[ev_id_it->first] = 0;
                ev_ids[ev_id_it->first] += ev_id_it->second;
            }

            for (ev_id_it = b->ev_ids.begin(); ev_id_it != b->ev_ids.end(); ++ev_id_it) {
                if ( ev_ids.find( ev_id_it->first ) == ev_ids.end() )
                    ev_ids[ev_id_it->first] = 0;
                ev_ids[ev_id_it->first] += ev_id_it->second;
            }

            weight += a->weight;
            weight += b->weight;
        }
    }
}
//}}}

//{{{ SV_BreakPoint:: ~SV_BreakPoint()
SV_BreakPoint::
~SV_BreakPoint()
{
    if (interval_l.p != NULL)
        free(interval_l.p);

    if (interval_r.p != NULL)
        free(interval_r.p);
}
//}}}

//{{{ void SV_BreakPoint:: free_evidence()
void
SV_BreakPoint::
free_evidence()
{
    vector<SV_Evidence*>::iterator it;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *sv_e = *it;
        delete sv_e;
    }
}
//}}}

//{{{ SV_BreakPoint:: SV_BreakPoint()
SV_BreakPoint::
SV_BreakPoint()
{
}
//}}}

//{{{ void SV_BreakPoint:: print_evidence()
void
SV_BreakPoint::
print_evidence(string pre)
{
    vector<SV_Evidence*>::iterator it;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *sv_e = *it;
        cout << pre;
        sv_e->print_evidence();
        // cout << pre;
        // SV_BreakPoint *tmp_bp = sv_e->get_bp();
        // tmp_bp->print_bedpe(-1,1);
        // delete(tmp_bp);
    }
}
//}}}

//{{{ ostream& operator << (ostream& out, const SV_BreakPoint& b)
ostream& operator<< (ostream& out, const SV_BreakPoint& b)
{
    out <<
        b.interval_l.i.chr << "," <<
        b.interval_l.i.start << "," <<
        (b.interval_l.i.end+1) << "," <<
        b.interval_l.i.strand << "\t" <<
        b.interval_r.i.chr << "," <<
        b.interval_r.i.start << "," <<
        (b.interval_r.i.end+1) << "," <<
        b.interval_r.i.strand << "\t" <<
        b.evidence.size();

    return out;
}
//}}}

//{{{ bool BreakPoint:: does_intersect(struct interval *a,
bool
SV_BreakPoint::
does_intersect(struct breakpoint_interval *a,
               struct breakpoint_interval *b,
               bool check_strand)
{
    if ( (a->i.chr == b->i.chr) &&
            (a->i.start < b->i.end) &&
            (a->i.end > b->i.start) )
        if (check_strand)
            return (a->i.strand == b->i.strand);
        else
            return true;
    else
        return false;
}
//}}}

//{{{ void SV_BreakPoint:: init_interval_probabilities()
void
SV_BreakPoint::
init_interval_probabilities()
{
    // skip if they have been initialized
    if ( (interval_l.p == NULL) || (interval_r.p == NULL) ) {

        // take the first piece of evidence from the list
        SV_Evidence *e = evidence[0];
        e->set_bp_interval_probability(&interval_l);
        e->set_bp_interval_probability(&interval_r);
    }
}
//}}}

//{{{ void SV_BreakPoint:: free_interval_probabilities()
void
SV_BreakPoint::
free_interval_probabilities()
{
    // skip if they have been initialized
    if ( (interval_l.p != NULL) || (interval_r.p != NULL) ) {

        // take the first piece of evidence from the list
        //SV_Evidence *e = evidence[0];
        free(interval_l.p);
        free(interval_r.p);
        interval_l.p = NULL;
        interval_r.p = NULL;
    }
}
//}}}

//{{{ void BreakPoint:: merge(BreakPoint *p)
bool
SV_BreakPoint::
merge(SV_BreakPoint *p)
{
    // p->a may not overlap a, so we need to check
    // after this a_overlap_intr will point to the interval in p that
    // overlaps the current a interval and bt_overlap_intr overlaps the
    // current b interval
    // if this->a overalps p->a then ll will be true, otherise false
    //bool ll = true;
    struct breakpoint_interval *a_overlap_intr, *b_overlap_intr;

    if ( does_intersect(&interval_l, &(p->interval_l), false) &&
            does_intersect(&interval_r, &(p->interval_r), false) ) {
        //ll = true;
        a_overlap_intr = &(p->interval_l);
        b_overlap_intr = &(p->interval_r);
    } else if ( does_intersect(&interval_r, &(p->interval_l), false) &&
                does_intersect(&interval_l, &(p->interval_r), false) ) {
        //ll = false;
        a_overlap_intr = &(p->interval_r);
        b_overlap_intr = &(p->interval_l);
    } else {
        /*
        cerr << "Error in merge(): no correct overlap" << endl <<
             "t:" <<
             interval_l.i.chr << "," <<
             interval_l.i.start << "," <<
             interval_l.i.end << "," <<
             interval_l.i.strand << " " <<
             interval_r.i.chr << "," <<
             interval_r.i.start << "," <<
             interval_r.i.end << "," <<
             interval_r.i.strand << "\t" <<
             "p:" <<
             p->interval_l.i.chr << "," <<
             p->interval_l.i.start << "," <<
             p->interval_l.i.end << "," <<
             p->interval_l.i.strand << " " <<
             p->interval_r.i.chr << "," <<
             p->interval_r.i.start << "," <<
             p->interval_r.i.end << "," <<
             p->interval_r.i.strand << "\t" <<
             "tl-pl:" <<
             does_intersect(&interval_l, &(p->interval_l), true) <<
             " " <<
             "tl-pr:" <<
             does_intersect(&interval_l, &(p->interval_r), true) <<
             " " <<
             "tr-pl:" <<
             does_intersect(&interval_r, &(p->interval_l), true) <<
             " " <<
             "tr-pr:" <<
             does_intersect(&interval_r, &(p->interval_r), true) <<
             endl;
        */
        return false;
    }

#if 1
    // just put everything together and refine the breakpoint boundaries to be
    // the mean of the evidence set
    vector<SV_Evidence*>::iterator ev_it;

    // copy all of the evidence and ids from the incomming breakpoint to the
    // current breakpoint
    for (ev_it = p->evidence.begin(); ev_it < p->evidence.end(); ++ev_it)
        evidence.push_back(*ev_it);

    map<int, int>::iterator ev_id_it;

    for (ev_id_it = p->ev_ids.begin(); ev_id_it != p->ev_ids.end(); ++ev_id_it) {
        if ( ev_ids.find( ev_id_it->first ) == ev_ids.end() )
            ev_ids[ev_id_it->first] = 0;

        ev_ids[ev_id_it->first] += ev_id_it->second;
    }

    // find the mean start and end points for the current bp
    CHR_POS l_start_sum = 0,
            l_end_sum = 0,
            r_start_sum = 0,
            r_end_sum = 0;
    for (ev_it = evidence.begin(); ev_it < evidence.end(); ++ev_it) {
        SV_Evidence *sv_e = *ev_it;
        SV_BreakPoint *tmp_bp = sv_e->get_bp();

        l_start_sum+=tmp_bp->interval_l.i.start;
        l_end_sum+=tmp_bp->interval_l.i.end;
        r_start_sum+=tmp_bp->interval_r.i.start;
        r_end_sum+=tmp_bp->interval_r.i.end;
        delete(tmp_bp);
    }

    CHR_POS l_merged_start, l_merged_end, r_merged_start, r_merged_end;


    l_merged_start = l_start_sum / evidence.size();
    l_merged_end = l_end_sum / evidence.size();
    r_merged_start = r_start_sum / evidence.size();
    r_merged_end = r_end_sum / evidence.size();

    /*
    cerr << "l_merged_start:" << l_merged_start << " " <<
            "l_merged_end:" << l_merged_end << "\t" <<
            "r_merged_start:" << r_merged_start << " " <<
            "r_merged_end:" << r_merged_end << endl;
    */

    interval_l.i.start = l_merged_start;
    interval_l.i.end = l_merged_end;

    interval_r.i.start = r_merged_start;
    interval_r.i.end = r_merged_end;

    this->weight += p->weight;

    return true;
#endif


#if 0
    // just put everything together and keep the first breakpoint position
//{{{
    vector<SV_Evidence*>::iterator ev_it;

    for (ev_it = p->evidence.begin(); ev_it < p->evidence.end(); ++ev_it)
        evidence.push_back(*ev_it);

    map<int, int>::iterator ev_id_it;

    for (ev_id_it = p->ev_ids.begin(); ev_id_it != p->ev_ids.end(); ++ev_id_it) {
        if ( ev_ids.find( ev_id_it->first ) == ev_ids.end() )
            ev_ids[ev_id_it->first] = 0;

        ev_ids[ev_id_it->first] += ev_id_it->second;
    }

    this->weight += p->weight;

    return true;
//}}}
#endif

#if 0
    // merge thow read and shrink the vector to be the rejoins common to both
    // end points
//{{{
    CHR_POS a_merged_start,
            a_merged_end,
            b_merged_start,
            b_merged_end;

    log_space *a_merged_prob, *b_merged_prob;

    bool a_merges = test_interval_merge(&interval_l,
                                        a_overlap_intr,
                                        &a_merged_start,
                                        &a_merged_end,
                                        &a_merged_prob);


    bool b_merges = test_interval_merge(&interval_r,
                                        b_overlap_intr,
                                        &b_merged_start,
                                        &b_merged_end,
                                        &b_merged_prob);

    if (a_merges && b_merges) {

#ifdef TRACE
        //{{{
        cout <<
             interval_l.i.start << ":" <<
             interval_l.i.end << " " <<
             ascii_interval_prob(&interval_l) << endl <<

             a_overlap_intr->i.start << ":" <<
             a_overlap_intr->i.end << " " <<
             ascii_interval_prob(a_overlap_intr) << endl<<


             a_merged_start << ":" <<
             a_merged_end << " " <<
             ascii_prob(a_merged_prob, a_merged_end - a_merged_start) << endl <<

             interval_r.i.start << ":"	 <<
             interval_r.i.end << " "	 <<

             b_overlap_intr->i.start << ":" <<
             b_overlap_intr->i.end << " -> " <<

             b_merged_start << ":" <<
             b_merged_end << " " <<

             ascii_interval_prob(&interval_r) << "\t" <<

             endl;
        //}}}
#endif

        interval_l.i.start = a_merged_start;
        interval_l.i.end = a_merged_end;
        free(interval_l.p);
        interval_l.p = a_merged_prob;

        interval_r.i.start = b_merged_start;
        interval_r.i.end = b_merged_end;
        free(interval_r.p);
        interval_r.p = b_merged_prob;

        vector<SV_Evidence*>::iterator ev_it;

        for (ev_it = p->evidence.begin(); ev_it < p->evidence.end(); ++ev_it)
            evidence.push_back(*ev_it);

        map<int, int>::iterator ev_id_it;

        for (ev_id_it = p->ev_ids.begin(); ev_id_it != p->ev_ids.end(); ++ev_id_it) {
            if ( ev_ids.find( ev_id_it->first ) == ev_ids.end() )
                ev_ids[ev_id_it->first] = 0;

            ev_ids[ev_id_it->first] += ev_id_it->second;
        }

        this->weight += p->weight;

        return true;
    } else {
        free(a_merged_prob);
        free(b_merged_prob);

        return false;
    }
//}}}
#endif
}
//}}}

//{{{bool test_interval_merge(struct interval *curr_intr,
bool
SV_BreakPoint::
test_interval_merge(struct breakpoint_interval *curr_intr,
                    struct breakpoint_interval *new_intr,
                    CHR_POS *merged_start,
                    CHR_POS *merged_end,
                    log_space **merged_prob)
{

    log_space ls_threshold = get_ls(p_merge_threshold);

    // Find how much to clip from the start of the current break point (this)
    // or the new break point (p)
    // Case 1:
    // curr:   |---------------...
    // new:        |-----------...
    // c_clip  ||||
    // Case 2:
    // curr:       |-------...
    // new:    |-----------...
    // n_clip  ||||


    CHR_POS curr_clip_start = 0;
    CHR_POS new_clip_start = 0;
    if (curr_intr->i.start < new_intr->i.start) // Case 1
        curr_clip_start = new_intr->i.start - curr_intr->i.start;
    else if (curr_intr->i.start > new_intr->i.start) // Case 2
        new_clip_start = curr_intr->i.start - new_intr->i.start;

    // Find how much to clip from the end of the current break point (this)
    // or the new bre ak point (p)
    // Case 1:
    // curr:   ...--------|
    // new:    ...----|
    // c_clip          ||||
    // Case 2:
    // curr:   ...----|
    // new:    ...--------|
    // n_clip          ||||
    CHR_POS curr_clip_end = 0;
    //CHR_POS new_clip_end = 0;
    if (curr_intr->i.end > new_intr->i.end) // Case 1
        curr_clip_end = curr_intr->i.end - new_intr->i.end;
    //else if (curr_intr->i.end < new_intr->i.end) // Case 2
    //new_clip_end = new_intr->i.end - curr_intr->i.end;


    // Length of the new break point interval
    unsigned int new_len = (curr_intr->i.end - curr_clip_end) -
                           (curr_intr->i.start + curr_clip_start) + 1;

    // Before we acutally merge these two breakpoints, we need to test if the
    // merged probability is greater than zero, if it is not, then we should
    // not merge the two
    unsigned int i = 0;

    *merged_start = curr_intr->i.start + curr_clip_start;
    *merged_end = curr_intr->i.end - curr_clip_end;

    *merged_prob = (log_space *) malloc(new_len * sizeof(log_space));

    log_space ls_max = -INFINITY;

    for (i = 0; i < new_len; ++i) {
        (*merged_prob)[i] = ls_multiply(curr_intr->p[i + curr_clip_start],
                                        new_intr->p[i + new_clip_start]);
    }


    for (i = 0; i < new_len; ++i)
        if ( (*merged_prob)[i] > ls_max )
            ls_max = (*merged_prob)[i];

    if (ls_max < ls_threshold )
        return false;
    else
        return true;

}
//}}}

//{{{ void BreakPoint::trim_intervals(double mean, double stdev, double
int
SV_BreakPoint::
trim_intervals()
{
    vector<SV_Evidence*>::iterator it;
    vector<SV_BreakPoint *> bps;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *e = *it;
        SV_BreakPoint *tmp_bp = e->get_bp();
        tmp_bp->init_interval_probabilities();

        bps.push_back(tmp_bp);
    }

#if 0
    CHR_POS start_l = UINT_MAX,
            start_r = UINT_MAX,
            end_l = 0,
            end_r = 0;

    log_space *l, *r;
    get_mixture(bps, &start_l, &start_r, &end_l, &end_r, &l, &r);
#endif

    CHR_POS p_start_l = UINT_MAX,
            p_start_r = UINT_MAX,
            p_end_l = 0,
            p_end_r = 0;
    log_space *p_l, *p_r;
    int pr = get_product(bps,
                        &p_start_l,
                        &p_start_r,
                        &p_end_l,
                        &p_end_r,
                        &p_l,
                        &p_r);
    if (pr == 0) {
        return 0;
    }

    CHR_POS p_l_size = p_end_l - p_start_l + 1;
    log_space *p_t = (log_space *) malloc(p_l_size * sizeof(log_space));
    normalize_ls(p_l_size, p_l, p_t);

    int p_l_trim_start, p_l_trim_end;
    trim_interval(p_t, p_l_size, &p_l_trim_start, &p_l_trim_end);


    interval_l.i.start = p_start_l + p_l_trim_start;
    interval_l.i.end = p_start_l + p_l_trim_end - 1;

    free(p_t);

    CHR_POS p_r_size = p_end_r - p_start_r + 1;
    p_t = (log_space *) malloc(p_r_size * sizeof(log_space));
    normalize_ls(p_r_size, p_r, p_t);
    int p_r_trim_start, p_r_trim_end;
    trim_interval(p_t, p_r_size, &p_r_trim_start, &p_r_trim_end);


    free(p_t);

    interval_r.i.start = p_start_r + p_r_trim_start;
    interval_r.i.end = p_start_r + p_r_trim_end - 1;


    free(p_l);
    free(p_r);

    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;
        delete tmp_bp;
    }


    return 1;
}
//}}}

//{{{ void BreakPoint:: trim_interval(struct interval *curr_interval,
void
SV_BreakPoint:: trim_interval(log_space *interval_v,
                              unsigned int size,
                              int *trim_start,
                              int *trim_end)
{
    log_space max = -INFINITY;
    unsigned int i, max_i;
    for (i = 0; i < size; ++i) {
        if (interval_v[i] > max) {
            max = interval_v[i];
            max_i = i;
        }
    }

    float sum = 0;
    for (i = 0; i < size; ++i)
        sum = sum + get_p(interval_v[i]);


    log_space total = max;
    unsigned int l = max_i, r = max_i;
    /*
     * starting at max_i, grow the [l,r] region to include the adjacnet
     * position with the largest value until the cummulative value is greater
     * than the given threshold
     */

    while ( ((l > 0) || (r < (size-1))) && (total < get_ls(p_trim_threshold)) ){
        if ( l == 0 ) {
            total = ls_add(total, interval_v[r+1]);
            ++r;
        } else if ( r == (size - 1) ) {
            total = ls_add(total, interval_v[l-1]);
            --l;
        } else if ( interval_v[l] > interval_v[r+1] ) {
            total = ls_add(total, interval_v[l-1]);
            --l;
        } else {
            total = ls_add(total,interval_v[r+1]);
            ++r;
        }
    }

    if ((l == r) && (l == max_i))
        r+=1;

    *trim_start = l;
    *trim_end = r;
}

//}}}

//{{{ string SV_BreakPoint:: ascii_interval_prob(struct breakpoint_interval *i)
string
SV_BreakPoint::
ascii_interval_prob(struct breakpoint_interval *i)
{
    int size = i->i.end - i->i.start;
    int j;
    stringstream oss;
    for (j = 0; j < size; ++j) {
        if (j!=0)
            oss << ",";
        oss << get_p(i->p[j]);
    }

    return oss.str();
}
//}}}

//{{{ string SV_BreakPoint:: ascii_prob(double *d, int size)
string
SV_BreakPoint::
ascii_prob(log_space *d, int size)
{
    int j;
    stringstream oss;
    for (j = 0; j < size; ++j) {
        oss << d[j] << " ";
    }

    return oss.str();
}
//}}}

//{{{ void SV_BreakPoint:: print_interval_probabilities()
void
SV_BreakPoint::
print_interval_probabilities()
{
    vector<SV_Evidence*>::iterator it;
    vector<SV_BreakPoint *> bps;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *e = *it;
        SV_BreakPoint *tmp_bp = e->get_bp();
        tmp_bp->init_interval_probabilities();

        bps.push_back(tmp_bp);
    }

    CHR_POS start_l = UINT_MAX,
            start_r = UINT_MAX,
            end_l = 0,
            end_r = 0;

    log_space *l, *r;
    //get_mixture(bps, &start_l, &start_r, &end_l, &end_r, &l, &r);
    get_product(bps, &start_l, &start_r, &end_l, &end_r, &l, &r);

    CHR_POS start_l_diff = interval_l.i.start - start_l,
            end_l_diff = end_l - interval_l.i.end,
            start_r_diff = interval_r.i.start - start_r,
            end_r_diff = end_r - interval_r.i.end;

    CHR_POS l_size = end_l - start_l + 1;
    log_space *t = (log_space *) malloc(l_size * sizeof(log_space));
    normalize_ls(l_size, l, t);

    CHR_POS i;

    cout << "PROBS:";

    for (i = start_l_diff; i < l_size - end_l_diff; i++) {
        if (i != start_l_diff)
            cout << ",";
        cout << get_p(t[i]);
    }

    cout << ";";

    free(t);

    CHR_POS r_size = end_r - start_r + 1;
    t = (log_space *) malloc(r_size * sizeof(log_space));
    normalize_ls(r_size, r, t);

    for (i = start_r_diff; i < r_size - end_r_diff; i++) {
        if (i != start_r_diff)
            cout << ",";
        cout << get_p(t[i]);
    }

    cout << "\t";

    free(t);


    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;
        delete tmp_bp;
    }

    free(l);
    free(r);
}
//}}}

//{{{ void SV_BreakPoint:: get_interval_probabilities(log_space *l,
void
SV_BreakPoint::
get_interval_probabilities(CHR_POS *start_l,
                           CHR_POS *start_r,
                           CHR_POS *end_l,
                           CHR_POS *end_r,
                           log_space **l,
                           log_space **r)
{
    vector<SV_Evidence*>::iterator it;
    vector<SV_BreakPoint *> bps;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *e = *it;
        SV_BreakPoint *tmp_bp = e->get_bp();
        tmp_bp->init_interval_probabilities();

        bps.push_back(tmp_bp);
    }

    /*
    CHR_POS start_l = UINT_MAX,
            start_r = UINT_MAX,
            end_l = 0,
            end_r = 0;
    */

    *start_l = UINT_MAX;
    *start_r = UINT_MAX;
    *end_l = 0;
    *end_r = 0;
    log_space *t_l, *t_r;
    //get_mixture(bps, &start_l, &start_r, &end_l, &end_r, &l, &r);
    //get_product(bps, &start_l, &start_r, &end_l, &end_r, &t_l, &t_r);
    get_product(bps, start_l, start_r, end_l, end_r, &t_l, &t_r);

    CHR_POS l_size = *end_l - *start_l + 1;
    *l = (log_space *) malloc(l_size * sizeof(log_space));
    normalize_ls(l_size, t_l, *l);

    free(t_l);

    CHR_POS r_size = *end_r - *start_r + 1;
    *r = (log_space *) malloc(r_size * sizeof(log_space));
    normalize_ls(r_size, t_r, *r);

    free(t_r);

    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;
        delete tmp_bp;
    }
}
//}}}

//{{{ void SV_BreakPoint:: print_bedpe(int id)
void
SV_BreakPoint::
print_bedpe(int id, int print_prob)
{
    map<string,int> uniq_strands;
    vector<SV_Evidence*>::iterator it;
    vector<SV_BreakPoint *> bps;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *e = *it;
        SV_BreakPoint *tmp_bp = e->get_bp();
        tmp_bp->init_interval_probabilities();

        bps.push_back(tmp_bp);

        stringstream strands;
        strands << tmp_bp->interval_l.i.strand << tmp_bp->interval_r.i.strand;

        if (uniq_strands.find(strands.str()) == uniq_strands.end())
            uniq_strands[strands.str()] = 1;
        else
            uniq_strands[strands.str()] = uniq_strands[strands.str()] + 1;
    }

    double score_l, score_r;
    get_score(bps, &score_l, &score_r);

    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;
        delete tmp_bp;
    }

    // use the address of the current object as the id
    string sep = "\t";

    CHR_POS open_l_start = 0,
            open_r_start = 0;

    if (interval_l.i.start > 0)
        open_l_start = interval_l.i.start - 1;

    if (interval_r.i.start > 0)
        open_r_start = interval_r.i.start - 1;


    cout <<
         interval_l.i.chr << sep <<
         open_l_start << sep <<
         interval_l.i.end  << sep <<
         interval_r.i.chr << sep <<
         open_r_start << sep<<
         interval_r.i.end << sep;

    if (id == -1)
        cout << this << sep;
    else
        cout << id << sep;

    cout <<
         (score_l+score_r) << "\t" <<
         interval_l.i.strand << "\t" <<
         interval_r.i.strand << "\t";

    cout << "TYPE:";

    if (interval_l.i.chr.compare(interval_r.i.chr) != 0)
        cout << "INTERCHROM";
    else if (type == DELETION )
        cout << "DELETION";
    else if (type == DUPLICATION)
        cout << "DUPLICATION";
    else if (type == INVERSION)
        cout << "INVERSION";
    else
        cout <<  "???";

    cout <<  "\t";

    map<int, int>::iterator ev_ids_it;
    vector<int> _ev_ids;
    for ( ev_ids_it = ev_ids.begin(); ev_ids_it != ev_ids.end(); ++ev_ids_it)
        _ev_ids.push_back(ev_ids_it->first);

    sort(_ev_ids.begin(), _ev_ids.end());

    vector<int>::iterator _ev_ids_it;
    map<string,int> ev_map;
    cout << "IDS:";
    for (_ev_ids_it = _ev_ids.begin();
	 _ev_ids_it != _ev_ids.end();
	 ++_ev_ids_it) {

	string samp = SV_EvidenceReader::sample_names[*_ev_ids_it];
    	string ev_type = SV_EvidenceReader::ev_types[*_ev_ids_it];

	ev_map[samp + "_" + ev_type] += ev_ids[*_ev_ids_it];
    }
    for (map<string,int>::iterator s_it = ev_map.begin();
	 s_it != ev_map.end();
	 ++s_it) {
        if (s_it != ev_map.begin())
            cout << ";";
        cout << s_it->first << "," << s_it->second;
    }

    cout << "\t";

    cout << "STRANDS:";
    map<string,int>:: iterator s_it;
    for ( s_it = uniq_strands.begin(); s_it != uniq_strands.end(); ++s_it) {
        if (s_it != uniq_strands.begin())
            cout <<  ";";
        cout << s_it->first << "," << s_it->second;
    }


    cout << "\t";


    // Get the most likley posistions

    CHR_POS start_l, start_r, end_l, end_r;
    log_space *l,*r;
    //get_interval_probabilities(&l,&r);

    get_interval_probabilities(&start_l,
                               &start_r,
                               &end_l,
                               &end_r,
                               &l,
                               &r);

    // l and r contain the full, untrimmed distribution.  All of the
    // calculations below must be shifted to the right by an offset so they are
    // considering the correction area of the pdf

    CHR_POS l_trim_offset = interval_l.i.start - start_l;
    CHR_POS r_trim_offset = interval_r.i.start - start_r;

    // find relative position of max value in left
    log_space max = -INFINITY;
    unsigned int i, l_max_i = 0, r_max_i = 0;
    for (i = 0; i < (interval_l.i.end - interval_l.i.start + 1); ++i) {
        if (l[l_trim_offset + i] > max) {
            max = l[l_trim_offset + i];
            l_max_i = l_trim_offset + i;
        }
    }

    // need to user start_l here, not interval_l.i.start, since l[0] is at
    // position start_l, and l_max_i is w.r.t. l[0]
    CHR_POS abs_max_l = start_l + l_max_i;

    // find relative position of max value in right
    max = -INFINITY;
    for (i = 0; i < (interval_r.i.end - interval_r.i.start + 1); ++i) {
        if (r[r_trim_offset + i] > max) {
            max = r[r_trim_offset + i];
            r_max_i = r_trim_offset + i;
        }
    }

    // need to user start_r here, not interval_r.i.start, since r[0] is at
    // position start_r, and r_max_i is w.r.t. r[0]
    CHR_POS abs_max_r = start_r + r_max_i;

    cout << "MAX:" << interval_l.i.chr << ":" << abs_max_l << ";"<<
        interval_r.i.chr << ":" << abs_max_r;

    cout << "\t";


    // Get the area that includes the max and 95% of the probabitliy
    log_space p_95 = get_ls(0.95);

    log_space total = l[l_max_i];
    CHR_POS l_l_i = l_max_i,
            l_r_i = l_max_i,
            t_last = interval_l.i.end - interval_l.i.start;

    while ( ((l_l_i > 0) || (l_r_i < t_last)) && (total < p_95) ){
        if ( l_l_i == 0 ) {
            total = ls_add(total, l[l_r_i+1]);;
            ++l_r_i;
        } else if ( l_r_i == t_last ) {
            total = ls_add(total, l[l_l_i-1]);
            --l_l_i;
        } else if ( l[l_l_i-1] == l[l_r_i+1] ) {
            total = ls_add(total, l[l_r_i+1]);
            total = ls_add(total, l[l_l_i-1]);
            --l_l_i;
            ++l_r_i;
        } else if ( l[l_l_i-1] > l[l_r_i+1] ) {
            total = ls_add(total, l[l_l_i-1]);
            --l_l_i;
        } else {
            total = ls_add(total,l[l_r_i+1]);
            ++l_r_i;
        }
    }


    CHR_POS abs_l_l_95 = start_l + l_l_i,
            abs_l_r_95 = start_l + l_r_i;

    total = r[r_max_i];
    CHR_POS r_l_i = r_max_i,
            r_r_i = r_max_i;
    t_last = interval_r.i.end - interval_r.i.start;

    while ( ((r_l_i > 0) || (r_r_i < t_last)) && (total < p_95) ){
        if ( r_l_i == 0 ) {
            total = ls_add(total, r[r_r_i+1]);
            ++r_r_i;
        } else if ( r_r_i == t_last ) {
            total = ls_add(total, r[r_l_i-1]);
            --r_l_i;
        } else if ( r[r_l_i-1] == r[r_r_i+1] ) {
            total = ls_add(total, r[r_r_i+1]);
            total = ls_add(total, r[r_l_i-1]);
            --r_l_i;
            ++r_r_i;
        } else if ( r[r_l_i-1] > r[r_r_i+1] ) {
            total = ls_add(total, r[r_l_i-1]);
            --r_l_i;
        } else {
            total = ls_add(total,r[r_r_i+1]);
            ++r_r_i;
        }
    }


    CHR_POS abs_r_l_95 = start_r + r_l_i,
            abs_r_r_95 = start_r + r_r_i;


    cout << "95:" <<
        interval_l.i.chr << ":" << abs_l_l_95 << "-"<< abs_l_r_95 <<";"<<
        interval_r.i.chr << ":" << abs_r_l_95 << "-"<< abs_r_r_95;

    if (print_prob > 0) {
        cout << "\tLP:";

        for (i = 0; i < (interval_l.i.end - interval_l.i.start + 1); ++i) {
            if (i != 0)
                cout << ",";
            cout << get_p(l[l_trim_offset+i]);
        }

        cout << "\tRP:";

        for (i = 0; i < (interval_r.i.end - interval_r.i.start + 1); ++i) {
            if (i != 0)
                cout << ",";
            cout << get_p(r[r_trim_offset+i]);
        }

    }

    free(l);
    free(r);
    cout << endl;
}
//}}}

//{{{void SV_BreakPoint:: cluster(UCSCBins<SV_BreakPoint*> &bins);
void
SV_BreakPoint::
cluster( UCSCBins<SV_BreakPoint*> &r_bin)
{
    vector< UCSCElement<SV_BreakPoint*> > tmp_hits_r =
        r_bin.get(interval_r.i.chr,
                  interval_r.i.start,
                  interval_r.i.end,
                  interval_r.i.strand,
                  false);

    if (tmp_hits_r.size() == 0) {
        insert(r_bin);
    } else {
        vector< UCSCElement<SV_BreakPoint*> >::iterator it;

        it = tmp_hits_r.begin();

        // remove hits that are of a different type, or do not intersect the
        // left side, hopefully there is only 1 guy left
        while(it != tmp_hits_r.end()) {
            if(it->value->type != type)
                it = tmp_hits_r.erase(it);
            //else if ( it->value->interval_l.i.chr != interval_l.i.chr )
            else if(it->value->interval_l.i.chr.compare(interval_l.i.chr) != 0)
                it = tmp_hits_r.erase(it);
            else if (!((it->value->interval_l.i.start < interval_l.i.end) &&
                       (it->value->interval_l.i.end > interval_l.i.start) ) )
                it = tmp_hits_r.erase(it);
            else
                ++it;
        }

        // non left to match against
        if (tmp_hits_r.size() < 1) {
            insert(r_bin);
        // one match so merge
        } else if  (tmp_hits_r.size() == 1) {
            // addr of the bp in both sets (and only item in the intersection)
            SV_BreakPoint *bp = tmp_hits_r[0].value;
            if ( this->merge(bp) ) {
                UCSCElement<SV_BreakPoint*> rm = tmp_hits_r[0];

                if (r_bin.remove(rm, false, false, true) != 0) {
                    cerr << "Error removing item in cluster()" << endl;
                    abort();
                }

                delete bp;
                this->insert(r_bin);
            } else {
                delete(this);
            }
            // more than one match
        } else {
            delete(this);
#if 0
//{{{
            cerr << "X" << endl;
            int top_hit = -1;
            int curr = 0;
            double top_score = -1;
            this->init_interval_probabilities();

            vector< UCSCElement<SV_BreakPoint*> >::iterator m_it;
            for(m_it = tmp_hits_r.begin(); m_it != tmp_hits_r.end(); ++m_it) {
                SV_BreakPoint *bp = m_it->value;

                struct breakpoint_interval *r_mixture, *l_mixture;

                bp->get_mixture(&l_mixture,&r_mixture);

                CHR_POS l_product_len, l_product_start, l_product_end;
                log_space *l_product_prob;

                interval_product(&interval_l,
                                 l_mixture,
                                 &l_product_len,
                                 &l_product_start,
                                 &l_product_end,
                                 &l_product_prob);

                log_space l_sum = get_ls(0);
                for (int i = 0; i < l_product_len; ++i)
                    l_sum = ls_add(l_sum, l_product_prob[i]);

                CHR_POS r_product_len, r_product_start, r_product_end;
                log_space *r_product_prob;

                interval_product(&interval_r,
                                 r_mixture,
                                 &r_product_len,
                                 &r_product_start,
                                 &r_product_end,
                                 &r_product_prob);

                log_space r_sum = get_ls(0);
                for (int i = 0; i < r_product_len; ++i)
                    r_sum = ls_add(r_sum, r_product_prob[i]);

                double score = (get_p(l_sum)+get_p(r_sum))/2;

                if (score > top_score) {
                    top_score = score;
                    top_hit = curr;
                }

                free(l_product_prob);

                free(r_mixture->p);
                free(l_mixture->p);
                free(r_mixture);
                free(l_mixture);

                ++curr;
            }

            this->free_interval_probabilities();

            SV_BreakPoint *bp = tmp_hits_r[top_hit].value;
            if ( this->merge(bp) ) {
                UCSCElement<SV_BreakPoint*> rm = tmp_hits_r[0];

                if (r_bin.remove(rm, false, false, true) != 0) {
                    cerr << "Error removing item in cluster()" << endl;
                    abort();
                }

                delete bp;
                this->insert(r_bin);
            }
//}}}
#endif
        }
    }
}
//}}}

//{{{ void SV_BreakPoint:: insert(UCSCBins<SV_BreakPoint*> &r_bin,
void
SV_BreakPoint::
insert(UCSCBins<SV_BreakPoint*> &r_bin)
{
    //cerr << "insert:"<< *this << endl;
    r_bin.add(interval_r.i.chr,
              interval_r.i.start,
              interval_r.i.end,
              interval_r.i.strand,
              this);
}
//}}}

//{{{ void SV_BreakPoint:: insert(UCSCBins<SV_BreakPoint*> &r_bin,
void
SV_BreakPoint::
insert(UCSCBins<SV_BreakPoint*> &l_bin,
       UCSCBins<SV_BreakPoint*> &r_bin)
{
    l_bin.add(interval_l.i.chr,
              interval_l.i.start,
              interval_l.i.end,
              interval_l.i.strand,
              this);

    r_bin.add(interval_r.i.chr,
              interval_r.i.start,
              interval_r.i.end,
              interval_r.i.strand,
              this);
}
//}}}

//{{{ vector<int> SV_BreakPoint:: get_evidence_ids()
vector<int>
SV_BreakPoint::
get_evidence_ids()
{
    vector<int> ev_ids;

    vector<SV_Evidence*>::iterator it;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *e = *it;
        int ev_id = e->ev_id;
        ev_ids.push_back(ev_id);
    }

    sort(ev_ids.begin(), ev_ids.end());

    vector<int>::iterator it_id;
    it_id = unique(ev_ids.begin(), ev_ids.end());

    ev_ids.resize( it_id - ev_ids.begin() );

    return ev_ids;
}
//}}}

//{{{ void SV_BreakPoint:: get_mixture()
void
SV_BreakPoint::
get_mixture(struct breakpoint_interval **l_mixture,
            struct breakpoint_interval **r_mixture)
{
    vector<SV_Evidence*>::iterator it;
    CHR_POS start_l = UINT_MAX,
            start_r = UINT_MAX,
            end_l = 0,
            end_r = 0;

    vector<SV_BreakPoint *> bps;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *e = *it;
        SV_BreakPoint *tmp_bp = e->get_bp();
        tmp_bp->init_interval_probabilities();

        bps.push_back(tmp_bp);
    }

    // find boundaries
    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;

        if (tmp_bp->interval_l.i.start < start_l)
            start_l = tmp_bp->interval_l.i.start;

        if (tmp_bp->interval_r.i.start < start_r)
            start_r = tmp_bp->interval_r.i.start;

        if (tmp_bp->interval_l.i.end > end_l)
            end_l = tmp_bp->interval_l.i.end;

        if (tmp_bp->interval_r.i.end > end_r)
            end_r = tmp_bp->interval_r.i.end;
    }

    // initialize the distribution
    log_space *l = (log_space *)malloc((end_l-start_l+1)*sizeof(log_space));
    for (CHR_POS i = 0; i <  end_l - start_l + 1; ++ i)
        l[i] = get_ls(0);

    log_space *r = (log_space *)malloc((end_r-start_r+1)*sizeof(log_space));
    for (CHR_POS i = 0; i <  end_r - start_r + 1; ++ i)
        r[i] = get_ls(0);

    // find mixture distribution
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;

        CHR_POS l_offset = tmp_bp->interval_l.i.start - start_l;
        CHR_POS l_size = tmp_bp->interval_l.i.end -
                         tmp_bp->interval_l.i.start + 1;

        log_space *t = (log_space *) malloc(l_size * sizeof(log_space));
        normalize_ls(l_size, tmp_bp->interval_l.p, t);

        for (CHR_POS i = 0; i < l_size; ++i)
            l[i+l_offset] = ls_add(l[i+l_offset], t[i]);

        free(t);

        CHR_POS r_offset = tmp_bp->interval_r.i.start - start_r;
        CHR_POS r_size = tmp_bp->interval_r.i.end -
                         tmp_bp->interval_r.i.start + 1;

        t = (log_space *) malloc(r_size * sizeof(log_space));
        normalize_ls(r_size, tmp_bp->interval_r.p, t);

        for (CHR_POS i = 0; i < r_size; ++i)
            r[i+r_offset] = ls_add(r[i+r_offset], t[i]);

        free(t);
    }

    log_space *t = (log_space *)malloc((end_l-start_l+1)*sizeof(log_space));
    normalize_ls(end_l-start_l+1, l, t);
    free(l);
    l = t;

    t = (log_space *)malloc((end_r-start_r+1)*sizeof(log_space));
    normalize_ls(end_r-start_r+1, r, t);
    free(r);
    r = t;

    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;
        delete tmp_bp;
    }

    *l_mixture = (struct breakpoint_interval *) malloc (
                     sizeof( struct breakpoint_interval));
    (*l_mixture)->i.start = start_l;
    (*l_mixture)->i.end = end_l;
    (*l_mixture)->p = l;

    *r_mixture = (struct breakpoint_interval *) malloc (
                     sizeof( struct breakpoint_interval));
    (*r_mixture)->i.start = start_r;
    (*r_mixture)->i.end = end_r;
    (*r_mixture)->p = r;
}
//}}}

//{{{ void SV_BreakPoint:: do_it()
void
SV_BreakPoint::
do_it()
{
    vector<SV_Evidence*>::iterator it;
    vector<SV_BreakPoint *> bps;
    for (it = evidence.begin(); it < evidence.end(); ++it) {
        SV_Evidence *e = *it;
        SV_BreakPoint *tmp_bp = e->get_bp();
        tmp_bp->init_interval_probabilities();

        bps.push_back(tmp_bp);
    }

    double score_l, score_r;
    get_score(bps, &score_l, &score_r);

    vector<SV_BreakPoint *>::iterator bp_it;

    int i;
    vector<int> to_rm;
    for (unsigned int i = 0; i < bps.size(); ++i) {
        double new_score_l, new_score_r;
        vector<SV_BreakPoint *> tmp_bps;
        for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
            tmp_bps.push_back(*bp_it);
        }
        tmp_bps.erase(tmp_bps.begin() + i);

        // if new_score_l > score_l then the cluster is better off without that
        // particualr piece of evidence
        get_score(tmp_bps, &new_score_l, &new_score_r);

        //cout << "\t" << (new_score_l-score_l) <<
        //"\t" << (new_score_r-score_r) << endl;

        if ( ((new_score_l-score_l) > 0 ) &&
                ((new_score_r-score_r) > 0 ) ) {
            to_rm.push_back(i);
            //cout << "\t-" << endl;
        }
    }

    for (i = to_rm.size(); i > 0; --i) {
        //cout << i << endl;
        bps.erase(bps.begin()+i);
        --weight;
    }

    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;
        delete tmp_bp;
    }

    //free(l);
    //free(r);
}
//}}}

//{{{void SV_BreakPoint:: get_mixture(vector<SV_BreakPoint *> &bps,
void
SV_BreakPoint::
get_mixture(vector<SV_BreakPoint *> &bps,
            CHR_POS *start_l,
            CHR_POS *start_r,
            CHR_POS *end_l,
            CHR_POS *end_r,
            log_space **l,
            log_space **r)
{
    *start_l = UINT_MAX;
    *start_r = UINT_MAX,
     *end_l = 0;
    *end_r = 0;

    // find boundaries
    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;

        if (tmp_bp->interval_l.i.start < *start_l)
            *start_l = tmp_bp->interval_l.i.start;

        if (tmp_bp->interval_r.i.start < *start_r)
            *start_r = tmp_bp->interval_r.i.start;

        if (tmp_bp->interval_l.i.end > *end_l)
            *end_l = tmp_bp->interval_l.i.end;

        if (tmp_bp->interval_r.i.end > *end_r)
            *end_r = tmp_bp->interval_r.i.end;
    }

    // initialize the distribution
    *l = (log_space *)malloc((*end_l - *start_l + 1)*sizeof(log_space));
    for (CHR_POS i = 0; i <  *end_l - *start_l + 1; ++ i)
        (*l)[i] = get_ls(0);

    *r = (log_space *)malloc((*end_r - *start_r + 1)*sizeof(log_space));
    for (CHR_POS i = 0; i <  *end_r - *start_r + 1; ++ i)
        (*r)[i] = get_ls(0);

    // find mixture distribution
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;

        CHR_POS l_offset = tmp_bp->interval_l.i.start - *start_l;
        CHR_POS l_size = tmp_bp->interval_l.i.end -
                         tmp_bp->interval_l.i.start + 1;

        log_space *t = (log_space *) malloc(l_size * sizeof(log_space));
        normalize_ls(l_size, tmp_bp->interval_l.p, t);

        for (CHR_POS i = 0; i < l_size; ++i)
            (*l)[i+l_offset] = ls_add((*l)[i+l_offset], t[i]);

        free(t);

        CHR_POS r_offset = tmp_bp->interval_r.i.start - *start_r;
        CHR_POS r_size = tmp_bp->interval_r.i.end -
                         tmp_bp->interval_r.i.start + 1;

        t = (log_space *) malloc(r_size * sizeof(log_space));
        normalize_ls(r_size, tmp_bp->interval_r.p, t);

        for (CHR_POS i = 0; i < r_size; ++i) {
            (*r)[i+r_offset] = ls_add((*r)[i+r_offset], t[i]);
        }

        free(t);
    }
}
//}}}

//{{{void SV_BreakPoint:: get_product(vector<SV_BreakPoint *> &bps,
int
SV_BreakPoint::
get_product(vector<SV_BreakPoint *> &bps,
            CHR_POS *start_l,
            CHR_POS *start_r,
            CHR_POS *end_l,
            CHR_POS *end_r,
            log_space **l,
            log_space **r)
{


    // check to see if all of the intervals intersect, if they do find the
    // product directly, if they don't then identify the region based on those
    // that intersect the maximum point from the mixture distribution

    CHR_POS t_end_l = UINT_MAX,
            t_end_r = UINT_MAX,
            t_start_l = 0,
            t_start_r = 0;

    // find the min common area
    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
        SV_BreakPoint *tmp_bp = *bp_it;

        // left side

        if (tmp_bp->interval_l.i.start > t_start_l)
            t_start_l = tmp_bp->interval_l.i.start;
        if (tmp_bp->interval_l.i.end < t_end_l)
            t_end_l = tmp_bp->interval_l.i.end;

        // right side
        if (tmp_bp->interval_r.i.start > t_start_r)
            t_start_r = tmp_bp->interval_r.i.start;
        if (tmp_bp->interval_r.i.end < t_end_r)
            t_end_r = tmp_bp->interval_r.i.end;
    }

    // some of the intervals did not intersect
    if ( (t_start_l >= t_end_l) || (t_start_r >= t_end_r) ) {

        CHR_POS m_start_l = UINT_MAX,
                m_start_r = UINT_MAX,
                m_end_l = 0,
                m_end_r = 0;

        log_space *m_l, *m_r;

        // get mixture with absolute possitions
        // [m_start_l,m_end_l] and [m_start_r, m_end_r]
        get_mixture(bps,
                    &m_start_l,
                    &m_start_r,
                    &m_end_l,
                    &m_end_r,
                    &m_l,
                    &m_r);

        // find relative position of max value in left
        log_space max = -INFINITY;
        unsigned int i, l_max_i, r_max_i;
        for (i = 0; i < (m_end_l - m_start_l); ++i) {
            if (m_l[i] > max) {
                max = m_l[i];
                l_max_i = i;
            }
        }
        CHR_POS abs_max_l = m_start_l + l_max_i;


        // find relative position of max value in right
        max = -INFINITY;
        for (i = 0; i < (m_end_r - m_start_r); ++i) {
            if (m_r[i] > max) {
                max = m_r[i];
                r_max_i = i;
            }
        }

        CHR_POS abs_max_r = m_start_r + r_max_i;

        // l_max_i is the relative poss of the max element in l
        // the absolute position of l_max_i is m_start_l + l_max_i

        // r_max_i is the relative poss of the max element in r
        // the absolute position of r_max_i is m_start_r + r_max_i

        *end_l = UINT_MAX;
        *end_r = UINT_MAX,
        *start_l = 0;
        *start_r = 0;

        int was_set = 0;

        // find the min common area
        //vector<SV_BreakPoint *>::iterator bp_it;
        for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
            SV_BreakPoint *tmp_bp = *bp_it;

            // make sure the current breakpoint intervals intersection both the
            // left max (m_start_l + l_max_i) and the right max
            // (m_start_r + r_max_i)
            if ( // left side intersect the abs pos of the max
                ((tmp_bp->interval_l.i.start <= abs_max_l) &&
                 (tmp_bp->interval_l.i.end >= abs_max_l)) &&
                // right side intersects the abs pos of max
                ((tmp_bp->interval_r.i.start <= abs_max_r) &&
                 (tmp_bp->interval_r.i.end >= abs_max_r)) ) {

                was_set = 1;

                // left side

                if (tmp_bp->interval_l.i.start > *start_l)
                    *start_l = tmp_bp->interval_l.i.start;
                if (tmp_bp->interval_l.i.end < *end_l)
                    *end_l = tmp_bp->interval_l.i.end;


                // right side
                if (tmp_bp->interval_r.i.start > *start_r)
                    *start_r = tmp_bp->interval_r.i.start;
                if (tmp_bp->interval_r.i.end < *end_r)
                    *end_r = tmp_bp->interval_r.i.end;
            }
        }

        if (was_set == 0)
            return 0;

        *l = (log_space *)malloc((*end_l - *start_l + 1)*sizeof(log_space));
        for (CHR_POS i = 0; i <  (*end_l - *start_l + 1); ++ i)
            (*l)[i] = get_ls(1);

        *r = (log_space *)malloc((*end_r - *start_r + 1)*sizeof(log_space));
        for (CHR_POS i = 0; i <  (*end_r - *start_r + 1); ++ i)
            (*r)[i] = get_ls(1);

        // find product
        for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
            SV_BreakPoint *tmp_bp = *bp_it;
            // make sure the current breakpoint intervals intersection both the
            // left max (m_start_l + l_max_i) and the right max
            // (m_start_r + r_max_i)
            if ( // left side intersect the abs pos of the max
                ((tmp_bp->interval_l.i.start <= abs_max_l) &&
                 (tmp_bp->interval_l.i.end >= abs_max_l)) &&
                // right side intersects the abs pos of max
                ((tmp_bp->interval_r.i.start <= abs_max_r) &&
                 (tmp_bp->interval_r.i.end >= abs_max_r)) ) {

                CHR_POS l_offset = *start_l - tmp_bp->interval_l.i.start;

                CHR_POS l_size = tmp_bp->interval_l.i.end -
                                 tmp_bp->interval_l.i.start + 1;
                log_space *t = (log_space *) malloc(l_size * sizeof(log_space));
                normalize_ls(l_size, tmp_bp->interval_l.p, t);

                for (CHR_POS i = 0; i <  (*end_l - *start_l + 1); ++ i)
                    (*l)[i] = ls_multiply((*l)[i], t[i + l_offset]);

                free(t);

                CHR_POS r_offset = *start_r - tmp_bp->interval_r.i.start;

                CHR_POS r_size = tmp_bp->interval_r.i.end -
                                 tmp_bp->interval_r.i.start + 1;
                t = (log_space *) malloc(r_size * sizeof(log_space));
                normalize_ls(r_size, tmp_bp->interval_r.p, t);

                for (CHR_POS i = 0; i <  (*end_r - *start_r + 1); ++ i)
                    (*r)[i] = ls_multiply((*r)[i], t[i + r_offset]);

                free(t);
            }
        }
    } else {

        *start_l = t_start_l;
        *end_l = t_end_l;

        *start_r = t_start_r;
        *end_r = t_end_r;

        *l = (log_space *)malloc((*end_l - *start_l + 1)*sizeof(log_space));
        for (CHR_POS i = 0; i <  (*end_l - *start_l + 1); ++ i)
            (*l)[i] = get_ls(1);

        *r = (log_space *)malloc((*end_r - *start_r + 1)*sizeof(log_space));
        for (CHR_POS i = 0; i <  (*end_r - *start_r + 1); ++ i)
            (*r)[i] = get_ls(1);

        // find product
        for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
            SV_BreakPoint *tmp_bp = *bp_it;

            CHR_POS l_offset = *start_l - tmp_bp->interval_l.i.start;

            CHR_POS l_size = tmp_bp->interval_l.i.end -
                             tmp_bp->interval_l.i.start + 1;
            log_space *t = (log_space *) malloc(l_size * sizeof(log_space));
            normalize_ls(l_size, tmp_bp->interval_l.p, t);

            for (CHR_POS i = 0; i <  (*end_l - *start_l + 1); ++ i)
                (*l)[i] = ls_multiply((*l)[i], t[i + l_offset]);

            free(t);

            CHR_POS r_offset = *start_r - tmp_bp->interval_r.i.start;

            CHR_POS r_size = tmp_bp->interval_r.i.end -
                             tmp_bp->interval_r.i.start + 1;
            t = (log_space *) malloc(r_size * sizeof(log_space));
            normalize_ls(r_size, tmp_bp->interval_r.p, t);

            for (CHR_POS i = 0; i <  (*end_r - *start_r + 1); ++ i)
                (*r)[i] = ls_multiply((*r)[i], t[i + r_offset]);

            free(t);
        }
    }

    return 1;
}
//}}}

// CC 2014-01-11: get_score is not used in the VCF
//{{{void SV_BreakPoint:: get_score( vector<SV_BreakPoint *> &bps,
void
SV_BreakPoint::
get_score( vector<SV_BreakPoint *> &bps,
           double *score_l,
           double *score_r)
{
    CHR_POS start_l = UINT_MAX,
            start_r = UINT_MAX,
            end_l = 0,
            end_r = 0;

    log_space *l, *r;
    //get_mixture(bps, &start_l, &start_r, &end_l, &end_r, &l, &r);
    get_product(bps, &start_l, &start_r, &end_l, &end_r, &l, &r);

    /*
    log_space sum_r = -INFINITY, sum_l = -INFINITY;

    for (CHR_POS i = 0; i <  end_l - start_l + 1; ++ i)
        sum_l = ls_add(sum_l, l[i]);

    for (CHR_POS i = 0; i <  end_r - start_r + 1; ++ i)
        sum_r = ls_add(sum_r, r[i]);
    */
    log_space max_r = -INFINITY, max_l = -INFINITY;

    for (CHR_POS i = 0; i <  end_l - start_l + 1; ++ i)
        if (l[i] > max_l)
            max_l = l[i];

    for (CHR_POS i = 0; i <  end_r - start_r + 1; ++ i)
        if (r[i] > max_r)
            max_r = l[i];


    double width_l = end_l - start_l + 1;
    double width_r = end_r - start_r + 1;


    *score_l = get_p(max_l)/width_l;
    *score_r = get_p(max_r)/width_r;
    free(l);
    free(r);
}
//}}}

//{{{pair<int,int> SV_BreakPoint:: min_pair( vector< vector<
pair<CHR_POS,CHR_POS>
SV_BreakPoint::
min_pair( vector< vector< pair<CHR_POS,CHR_POS> > > &m)
{
    CHR_POS min_dist = UINT_MAX;
    CHR_POS max_overlap = 0;
    CHR_POS min_row = 0, min_col = 0;
    CHR_POS curr_dist, curr_overlap;

    unsigned int r, c;
    for (r = 0; r < m.size(); ++r) {
        vector< pair<CHR_POS,CHR_POS> > row = m[r];
        for (c = r + 1; c < row.size(); ++c) {
            curr_dist = row[c].first;
            curr_overlap = row[c].second;
            if (min_dist > curr_dist) {
                min_dist = curr_dist;
                max_overlap = curr_overlap;
                min_row = r;
                min_col = c;
            } else if ( (min_dist == curr_dist) && (max_overlap < curr_overlap)) {
                min_dist = curr_dist;
                max_overlap = curr_overlap;
                min_row = r;
                min_col = c;
            }
        }
    }

    pair<CHR_POS,CHR_POS> min_index;
    min_index.first = min_row;
    min_index.second = min_col;

    return min_index;
}
//}}}

//{{{void SV_BreakPoint:: get_max_range(struct breakpoint_interval *b,
void
SV_BreakPoint::
get_max_range(struct breakpoint_interval *b,
              CHR_POS *max_start,
              CHR_POS *max_end,
              log_space *max_value)
{
    // find the range of values that are in the max position for each bp
    CHR_POS i;
    *max_start = 0;
    *max_end = 0;
    *max_value = -INFINITY;
    for (i = 0; i < b->i.end - b->i.start; ++i) {
        if ( b->p[i] > *max_value ) {
            *max_value = b->p[i];
            *max_start = i;
        } else if ( b->p[i] == *max_value ) {
            *max_end = i;
        }
    }

    if (*max_end < *max_start)
        *max_end = *max_start;

    *max_end = b->i.start + *max_end;
    *max_start = b->i.start + *max_start;
}
//}}}

//{{{bool SV_BreakPoint:: peak_distance(struct breakpoint_interval
void
SV_BreakPoint::
peak_distance(SV_BreakPoint *e,
              CHR_POS *dist,
              CHR_POS *overlap)
{
    init_interval_probabilities();
    e->init_interval_probabilities();

    CHR_POS a_max_start, a_max_end, b_max_start, b_max_end;
    log_space a_max, b_max;
    get_max_range(&interval_l, &a_max_start, &a_max_end, &a_max);
    get_max_range(&(e->interval_l), &b_max_start, &b_max_end, &b_max);

    free_interval_probabilities();
    e->free_interval_probabilities();

    if ( (a_max_end >= b_max_start) && (a_max_start <= b_max_end) ) {
        *dist = 0;
        *overlap = min(a_max_end,b_max_end) - max(a_max_start,b_max_start) + 1;
    } else if (a_max_start > b_max_end) {
        *dist = a_max_start - b_max_end;
        *overlap = 0;
    } else {
        *dist = b_max_start - a_max_end;
        *overlap = 0;
    }
#if 0
    log_space ls_threshold = get_ls(p_merge_threshold);

    // Find how much to clip from the start of the current break point (this)
    // or the new break point (p)
    // Case 1:
    // curr:   |---------------...
    // new:        |-----------...
    // c_clip  ||||
    // Case 2:
    // curr:       |-------...
    // new:    |-----------...
    // n_clip  ||||


    CHR_POS curr_clip_start = 0;
    CHR_POS new_clip_start = 0;
    if (curr_intr->i.start < new_intr->i.start) // Case 1
        curr_clip_start = new_intr->i.start - curr_intr->i.start;
    else if (curr_intr->i.start > new_intr->i.start) // Case 2
        new_clip_start = curr_intr->i.start - new_intr->i.start;

    // Find how much to clip from the end of the current break point (this)
    // or the new bre ak point (p)
    // Case 1:
    // curr:   ...--------|
    // new:    ...----|
    // c_clip          ||||
    // Case 2:
    // curr:   ...----|
    // new:    ...--------|
    // n_clip          ||||
    CHR_POS curr_clip_end = 0;
    //CHR_POS new_clip_end = 0;
    if (curr_intr->i.end > new_intr->i.end) // Case 1
        curr_clip_end = curr_intr->i.end - new_intr->i.end;
    //else if (curr_intr->i.end < new_intr->i.end) // Case 2
    //new_clip_end = new_intr->i.end - curr_intr->i.end;


    // Length of the new break point interval
    unsigned int new_len = (curr_intr->i.end - curr_clip_end) -
                           (curr_intr->i.start + curr_clip_start) + 1;

    // Before we acutally merge these two breakpoints, we need to test if the
    // merged probability is greater than zero, if it is not, then we should
    // not merge the two
    unsigned int i = 0;

    *merged_start = curr_intr->i.start + curr_clip_start;
    *merged_end = curr_intr->i.end - curr_clip_end;

    *merged_prob = (log_space *) malloc(new_len * sizeof(log_space));

    log_space ls_max = -INFINITY;

    for (i = 0; i < new_len; ++i) {
        (*merged_prob)[i] = ls_multiply(curr_intr->p[i + curr_clip_start],
                                        new_intr->p[i + new_clip_start]);
    }


    for (i = 0; i < new_len; ++i)
        if ( (*merged_prob)[i] > ls_max )
            ls_max = (*merged_prob)[i];

    if (ls_max < ls_threshold )
        return false;
    else
        return true;

#endif
}
//}}}

int
SV_BreakPoint::
get_max_sample_weight()
{
    map<string,int> sample_weights;
    int max_sample_weight = 0;

    map<int, int>::iterator ids_it;
    for (ids_it = this->ev_ids.begin();
	 ids_it != this->ev_ids.end();
	 ++ids_it)
	sample_weights[SV_EvidenceReader::sample_names[ids_it->first]] += ids_it->second;

    map<string,int>::iterator sw_it;
    for (sw_it = sample_weights.begin();
    	 sw_it != sample_weights.end();
    	 ++sw_it) {
    	if (sw_it->second > max_sample_weight)
    	    max_sample_weight = sw_it->second;
    }
    return max_sample_weight;
}

#if 0
//{{{ boost::numeric::ublas::matrix<double>* BreakPoint::get_matrix(double
//  number of rows equals the lenght of the a interval
//  number of cols equals the lenght of the b interval
//  m(a,b)
boost::numeric::ublas::matrix<log_space>*
SV_BreakPoint::
get_matrix()
{
    // initialize the matrix
    boost::numeric::ublas::matrix<log_space> *m  =
        new boost::numeric::ublas::matrix<log_space>(
        interval_l.i.end - interval_l.i.start + 1,
        interval_r.i.end - interval_r.i.start + 1);

    // Print the interval vectors that will be used to set the martix priors
    // set the matrix priors to be the probability based on the left and right
    // intervals
    for (unsigned int a = 0; a < m->size1 (); ++a)
        for (unsigned int b = 0; b < m->size2 (); ++b)
            (*m)(a, b) = ls_multiply(interval_l.p[a], interval_r.p[b]);

    // Print the matrix priors
#if 0
    cerr << m->size1() << "," << m->size2() << endl;
    for (unsigned int l = 0; l < m->size1 (); ++l)
        for (unsigned int r = 0; r < m->size2 (); ++r)
            cout << (*m)(l,r) << endl;
#endif

    // For each pair in this breakpoint, we need to update the matrix
    // probabilities based on the template length of the pair given
    // breakpoint
    // Use log-space probabilities to avoid underflow
    vector<SV_Evidence*>::iterator it;

    for (it = evidence.begin(); it < evidence.end(); ++it)
        (*it)->update_matrix(m);

    log_space ls_sum = -INFINITY;

    for (unsigned int a = 0; a < m->size1 (); ++ a)
        for (unsigned int b = 0; b < m->size2 (); ++ b)
            ls_sum = ls_add(ls_sum, (*m)(a, b));


    for (unsigned int a = 0; a < m->size1 (); ++ a)
        for (unsigned int b = 0; b < m->size2 (); ++ b)
            (*m)(a, b) = ls_divide( (*m)(a, b), ls_sum );

    // Print the matrix
#if 0
    cerr << m->size1() << "," << m->size2() << endl;
    for (unsigned int l = 0; l < m->size1 (); ++l)
        for (unsigned int r = 0; r < m->size2 (); ++r)
            cout << (*m)(l,r) << endl;
#endif
    return m;
}
///}}}

//{{{ void SV_BreakPoint:: get_distances(vector< SV_BreakPoint*> &new_v)
void
SV_BreakPoint::
get_distances(vector< SV_BreakPoint*> &new_v)
{
    vector<SV_BreakPoint*>::iterator it;
    for (it = new_v.begin(); it != new_v.end(); ++it) {
        SV_BreakPoint *b = *it;
        b->init_interval_probabilities();
    }

    for (it = new_v.begin(); it != new_v.end(); ++it) {
        vector<SV_BreakPoint*>::iterator jt;
        SV_BreakPoint *b_i = *it;
        ostringstream oss;
        for (jt = new_v.begin(); jt != new_v.end(); ++jt) {
            CHR_POS product_len, product_start, product_end;
            log_space *product_prob;
            SV_BreakPoint *b_j = *jt;
            interval_product(&(b_i->interval_l),
                             &(b_j->interval_l),
                             &product_len,
                             &product_start,
                             &product_end,
                             &product_prob);
            if (product_len > 0) {
                log_space distance = -INFINITY;
                for (unsigned int i = 0; i < product_len; ++i)
                    distance = ls_add(distance,product_prob[i]);
                oss << get_p(distance) << "\t";
                free(product_prob);
            }
        }
        cerr << oss.str() << endl;
    }

    for (it = new_v.begin(); it != new_v.end(); ++it) {
        SV_BreakPoint *b = *it;
        b->free_interval_probabilities();
    }


}
//}}}

//{{{void interval_product(struct interval *curr_intr,
void
SV_BreakPoint::
interval_product(struct breakpoint_interval *a_intr,
                 struct breakpoint_interval *b_intr,
                 CHR_POS *product_len,
                 CHR_POS *product_start,
                 CHR_POS *product_end,
                 log_space **product_prob)
{

    // if the those intervals do not intersect, then set the len to zero and
    // return
    if ( (a_intr->i.end < b_intr->i.start) ||
            (a_intr->i.start > b_intr->i.end) ) {
        *product_len = 0;
        return;
    }


    // Find how much to clip from the start of the current break point (this)
    // or the new break point (p)
    // Case 1:
    // curr:   |---------------...
    // new:        |-----------...
    // c_clip  ||||
    // Case 2:
    // curr:       |-------...
    // new:    |-----------...
    // n_clip  ||||


    CHR_POS curr_clip_start = 0;
    CHR_POS new_clip_start = 0;
    if (a_intr->i.start < b_intr->i.start) // Case 1
        curr_clip_start = b_intr->i.start - a_intr->i.start;
    else if (a_intr->i.start > b_intr->i.start) // Case 2
        new_clip_start = a_intr->i.start - b_intr->i.start;

    // Find how much to clip from the end of the current break point (this)
    // or the new bre ak point (p)
    // Case 1:
    // curr:   ...--------|
    // new:    ...----|
    // c_clip          ||||
    // Case 2:
    // curr:   ...----|
    // new:    ...--------|
    // n_clip          ||||
    CHR_POS curr_clip_end = 0;
    if (a_intr->i.end > b_intr->i.end) // Case 1
        curr_clip_end = a_intr->i.end - b_intr->i.end;


    // Length of the new break point interval
    unsigned int new_len = (a_intr->i.end - curr_clip_end) -
                           (a_intr->i.start + curr_clip_start) + 1;

    // Before we acutally merge these two breakpoints, we need to test if the
    // merged probability is greater than zero, if it is not, then we should
    // not merge the two
    unsigned int i = 0;

    *product_start = a_intr->i.start + curr_clip_start;
    *product_end = a_intr->i.end - curr_clip_end;

    log_space *curr_p = (log_space *)
                        malloc( (a_intr->i.end - a_intr->i.start + 1) *
                                sizeof(log_space));
    normalize_ls((a_intr->i.end - a_intr->i.start + 1),
                 a_intr->p,
                 curr_p);

    log_space *new_p = (log_space *)
                       malloc( (b_intr->i.end - b_intr->i.start + 1) *
                               sizeof(log_space));
    normalize_ls((b_intr->i.end - b_intr->i.start + 1),
                 b_intr->p,
                 new_p);

    *product_prob = (log_space *) malloc(new_len * sizeof(log_space));

    for (i = 0; i < new_len; ++i)
        (*product_prob)[i] = ls_multiply(a_intr->p[i + curr_clip_start],
                                         b_intr->p[i + new_clip_start]);

    *product_len = new_len;

    free(curr_p);
    free(new_p);
}
//}}}
#endif
