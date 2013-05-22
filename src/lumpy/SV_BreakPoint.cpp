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

	ids[e->id] = 1;
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

			map<int, int>::iterator id_it;

			for (id_it = a->ids.begin(); id_it != a->ids.end(); ++id_it) {
				if ( ids.find( id_it->first ) == ids.end() )
					ids[id_it->first] = 0;
				ids[id_it->first] += id_it->second;
			}

			for (id_it = a->ids.begin(); id_it != a->ids.end(); ++id_it) {
				if ( ids.find( id_it->first ) == ids.end() )
					ids[id_it->first] = 0;
				ids[id_it->first] += id_it->second;
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
		SV_BreakPoint *tmp_bp = sv_e->get_bp();
		tmp_bp->print_bedpe(-1);
		delete(tmp_bp);
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
		SV_Evidence *e = evidence[0];
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
		//{{{ Give error and abort
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
		//abort();
		return false;
		//}}}
	}

#if 1
	// just put everything together and refine the breakpoint boundaries to be
	// the mean of the evidence set
	vector<SV_Evidence*>::iterator ev_it;

	// copy all of the evidence and ids from the incomming breakpoint to the
	// current breakpoint 
	for (ev_it = p->evidence.begin(); ev_it < p->evidence.end(); ++ev_it)
		evidence.push_back(*ev_it);

	map<int, int>::iterator id_it;

	for (id_it = p->ids.begin(); id_it != p->ids.end(); ++id_it) {
		if ( ids.find( id_it->first ) == ids.end() )
			ids[id_it->first] = 0;

		ids[id_it->first] += id_it->second;
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

	map<int, int>::iterator id_it;

	for (id_it = p->ids.begin(); id_it != p->ids.end(); ++id_it) {
		if ( ids.find( id_it->first ) == ids.end() )
			ids[id_it->first] = 0;

		ids[id_it->first] += id_it->second;
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

		map<int, int>::iterator id_it;

		for (id_it = p->ids.begin(); id_it != p->ids.end(); ++id_it) {
			if ( ids.find( id_it->first ) == ids.end() )
				ids[id_it->first] = 0;

			ids[id_it->first] += id_it->second;
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
void
SV_BreakPoint::
trim_intervals()
{
	int a_v_size = interval_l.i.end - interval_l.i.start + 1,
		b_v_size = interval_r.i.end - interval_r.i.start + 1;
	log_space *a_v = (log_space *) malloc( a_v_size* sizeof(log_space));
	log_space *b_v = (log_space *) malloc( b_v_size* sizeof(log_space));

	for (int a = 0; a < a_v_size; ++a) 
		a_v[a] = interval_l.p[a];

	for (int b = 0; b < b_v_size; ++b)
		b_v[b] = interval_r.p[b];

	for (int b = 0; b < b_v_size; ++b)
		b_v[b] = interval_r.p[b];

	//cerr << ascii_interval_prob(&interval_l) <<endl;
	trim_interval(&interval_l, a_v, a_v_size);
	//cerr << ascii_interval_prob(&interval_l) <<endl;
	free(a_v);

	//cerr << ascii_interval_prob(&interval_r) <<endl;
	trim_interval(&interval_r, b_v, b_v_size);
	//cerr << ascii_interval_prob(&interval_r) <<endl;
	free(b_v);
}
//}}}

//{{{ void BreakPoint:: trim_interval(struct interval *curr_interval,
void
SV_BreakPoint:: trim_interval(struct breakpoint_interval *curr_interval,
							  log_space *interval_v,
							  unsigned int size)
{
	log_space max = -INFINITY;
	unsigned int i;
	for (i = 0; i < size; ++i)
		if (interval_v[i] > max)
			max = interval_v[i];
	
	int v_first = -1;
	for (i = 0; i < size; ++i) 
		if (get_p(interval_v[i])/get_p(max) < p_trim_threshold) 
			v_first=i;
		else
			break;

	int v_last = -1;
	for (i = size - 1; i > 0; --i)
		if (get_p(interval_v[i])/get_p(max) < p_trim_threshold) 
			v_last=i;
		else
			break;


	if ( (v_first != -1) && (v_last != -1) ) {
		// Adjust the breakpoint intervals
		curr_interval->i.start = curr_interval->i.start + v_first;
		curr_interval->i.end = curr_interval->i.start + (v_last - v_first);

		unsigned int v_size = curr_interval->i.end - curr_interval->i.start + 1;

		free(curr_interval->p);
		curr_interval->p = (log_space *) malloc(v_size * sizeof(log_space));

		for (i = 0; i < v_size; ++i) 
			curr_interval->p[i] = interval_v[i + v_first];
	}
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

//{{{ void SV_BreakPoint:: print_bedpe(int score)
void
SV_BreakPoint::
print_bedpe(int score)
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
	for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
		SV_BreakPoint *tmp_bp = *bp_it;
        delete tmp_bp;
    }


	// use the address of the current object as the id
	string sep = "\t";
	cout << 
		interval_l.i.chr << sep <<
		interval_l.i.start << sep <<
		(interval_l.i.end + 1) << sep <<
		interval_r.i.chr << sep <<
		interval_r.i.start << sep<<
		(interval_r.i.end + 1) << sep<<
		this << sep <<
		(score_l+score_r) << "\t" <<
		interval_l.i.strand << "\t" <<
		interval_r.i.strand << "\t";

		cout << "TYPE:";

		if (type == DELETION )
        	cout << "DELETION";
		else if (type == DUPLICATION)
        	cout << "DUPLICATION";
        else if (type == INVERSION)
        	cout << "INVERSION";
        else if (type == TRANSLOCATION)
        	cout <<  "TRANSLOCATION";
        else
        	cout <<  "???";

		cout <<  "\t";

		//ascii_interval_prob(&interval_l) << "\t" <<
		//ascii_interval_prob(&interval_r) <<

	/*
	vector <int> ids = get_evidence_ids();
	vector <int>::iterator i;
	cout << "ids:";
	for (i=ids.begin(); i!=ids.end(); ++i) {
		if (i != ids.begin())
			cout << ",";
		cout << *i;
	}
	*/

	map<int, int>::iterator ids_it;
	vector<int> _ids;
	for ( ids_it = ids.begin(); ids_it != ids.end(); ++ids_it)
		_ids.push_back(ids_it->first);

	sort(_ids.begin(), _ids.end());

	vector<int>::iterator _ids_it;

	cout << "IDS:";
	for ( _ids_it = _ids.begin(); _ids_it != _ids.end(); ++_ids_it) {
		if (_ids_it != _ids.begin())
			cout << ";";
		cout << *_ids_it << "," << ids[*_ids_it];
	}
		
	cout << endl;
}
//}}}

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
                for (int i = 0; i < product_len; ++i)
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
		if (tmp_hits_r.size() < 1)
			insert(r_bin);
		// one match so merge
		else if  (tmp_hits_r.size() == 1) {
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
			}
        // more than one match
		} else {
        #if 0
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
	//cerr << *this << endl;
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
	vector<int> ids;

	vector<SV_Evidence*>::iterator it;
	for (it = evidence.begin(); it < evidence.end(); ++it) {
		SV_Evidence *e = *it;
		int id = e->id;
		ids.push_back(id);
	}

	sort(ids.begin(), ids.end());

	vector<int>::iterator it_id;
	it_id = unique(ids.begin(), ids.end());

	ids.resize( it_id - ids.begin() );

	return ids;
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
    for (i = 0; i < bps.size(); ++i) {
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

        for (CHR_POS i = 0; i < r_size; ++i) {
            r[i+r_offset] = ls_add(r[i+r_offset], t[i]);
        }

        free(t);
    }

    log_space sum_r = -INFINITY, sum_l = -INFINITY;

    for (CHR_POS i = 0; i <  end_l - start_l + 1; ++ i)
        sum_l = ls_add(sum_l, l[i]);

    for (CHR_POS i = 0; i <  end_r - start_r + 1; ++ i)
        sum_r = ls_add(sum_r, r[i]);

    double width_l = end_l - start_l + 1;
    double width_r = end_r - start_r + 1;

    /*
    cerr << "\t\t" << max_l << "," << get_p(max_l) << "," << width_l << "\t" <<
        max_r << "," << get_p(max_r) << "," << width_r << endl;
    */

    *score_l = get_p(sum_l)/width_l;
    *score_r = get_p(sum_r)/width_r;
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

	int r, c;
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
			} else if ( (min_dist == curr_dist) && (max_overlap < curr_overlap)){
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
#endif

