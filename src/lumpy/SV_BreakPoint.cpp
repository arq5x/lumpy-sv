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

//{{{ SV_BreakPoint:: ~SV_BreakPoint()
SV_BreakPoint::
~SV_BreakPoint()
{
	free(interval_l.p);
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
		//sv_e->print_evidence();
		SV_BreakPoint *tmp_bp = sv_e->get_bp();
		tmp_bp->print_bedpe(-1);
		//tmp_bp->init_interval_probabilities();
		//struct breakpoint_interval i = tmp_bp->interval_l;
		//cout <<  ascii_interval_prob(&i) << endl;
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

//{{{ void BreakPoint:: merge(BreakPoint *p)
bool
SV_BreakPoint::
merge(SV_BreakPoint *p)
{
	//at this point malloc that space for the arrays in both this and in p
	init_interval_probabilities();
	p->init_interval_probabilities();


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

//{{{ void BreakPoint::trim_intervals(double mean, double stdev, double
void
SV_BreakPoint::
trim_intervals()
{
#if 0
	//cerr << "Get Matrix..." << endl;
	boost::numeric::ublas::matrix<log_space> *m = get_matrix();
	//cerr << "Done." << endl;

	//  Get the sum of probabilities for each row/col to get the probability
	//  for the left/right end
	 
	int a_v_size = m->size1(),
		b_v_size = m->size2();

	log_space *a_v = (log_space *) malloc(a_v_size * sizeof(log_space));
	log_space *b_v = (log_space *) malloc(b_v_size * sizeof(log_space));

	unsigned int a, b;

	for (a = 0; a < a_v_size; ++a) 
		a_v[a] = -INFINITY;

	for (b = 0; b < b_v_size; ++b)
		b_v[b] = -INFINITY;

	//cerr<< "Get Vectors..." << endl;
	for (unsigned int a = 0; a < a_v_size; ++a)
		for (unsigned int b = 0; b < b_v_size; ++b)
			a_v[a] = ls_add(a_v[a], (*m)(a,b));

	for (unsigned int a = 0; a < a_v_size; ++a)
		for (unsigned int b = 0; b < b_v_size; ++b)
			b_v[b] = ls_add(b_v[b], (*m)(a,b));
	//cerr << "Done" << endl;
#endif

	int a_v_size = interval_l.i.end - interval_l.i.start + 1,
		b_v_size = interval_r.i.end - interval_r.i.start + 1;
	log_space *a_v = (log_space *) malloc( a_v_size* sizeof(log_space));
	log_space *b_v = (log_space *) malloc( b_v_size* sizeof(log_space));

	for (int a = 0; a < a_v_size; ++a) 
		a_v[a] = interval_l.p[a];

	for (int b = 0; b < b_v_size; ++b)
		b_v[b] = interval_r.p[b];

	trim_interval(&interval_l, a_v, a_v_size);
	free(a_v);

	trim_interval(&interval_r, b_v, b_v_size);
	free(b_v);
}
//}}}

//{{{ void BreakPoint:: trim_interval(struct interval *curr_interval,
void
SV_BreakPoint:: trim_interval(struct breakpoint_interval *curr_interval,
							  log_space *interval_v,
							  unsigned int size)
{
	//log_space ls_t = get_ls(p_trim_threshold);

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
	for (i = size - 1; i >= 0; --i)
		if (get_p(interval_v[i])/get_p(max) < p_trim_threshold) 
			v_last=i;
		else
			break;

//	v_first = -1;
//	v_last = -1;
//	// Find the range of values that are above the threshold
//	for (i = 0; i < size; ++i) {
//		if (interval_v[i] > ls_t)
//			if (v_first > -1)
//				v_last = i;
//			else {
//				v_first = i;
//				v_last = i;
//			}
//	}

	//cerr << v_first << "," << v_last << endl;

	if ( (v_first != -1) && (v_last != -1) ) {

		//cerr << v_first << "\t" << v_last << endl;

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
		oss << i->p[j] << " ";
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
		weight << "\t" <<
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

//{{{void SV_BreakPoint:: cluster(UCSCBins<SV_BreakPoint*> &bins);
void
SV_BreakPoint::
cluster( UCSCBins<SV_BreakPoint*> &r_bin)
{
#if 1
	vector< UCSCElement<SV_BreakPoint*> > tmp_hits_r =
			r_bin.get(interval_r.i.chr,
					  interval_r.i.start,
					  interval_r.i.end,
					  interval_r.i.strand,
					  false);		

	if (tmp_hits_r.size() == 0) {
		//insert(l_bin,r_bin);
		insert(r_bin);
	} else { 

		vector< UCSCElement<SV_BreakPoint*> >::iterator it;
		it = tmp_hits_r.begin();

		// remove hits that are of a different type, or do not intersect the
		// left side, hopefully there is only 1 guy left
		while(it != tmp_hits_r.end()) {
			if(it->value->type != type)
				it = tmp_hits_r.erase(it);
			else if ( it->value->interval_l.i.chr != interval_l.i.chr )
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
		//	insert(l_bin,r_bin);
		// one match so merge
		else if  (tmp_hits_r.size() == 1) {
			// addr of the bp in both sets (and only item in the intersection)

			SV_BreakPoint *bp = tmp_hits_r[0].value;

			if ( bp->merge(this) ) {
				UCSCElement<SV_BreakPoint*> rm = tmp_hits_r[0];
				r_bin.remove(rm, false, false, true);

				delete this;

				//bp->insert(l_bin,r_bin);
				bp->insert(r_bin);

			}
		} else {
		}
	}
#endif

#if 0
	// Find matching break points
	vector< UCSCElement<SV_BreakPoint*> > tmp_hits_l =
			l_bin.get(interval_l.i.chr,
					  interval_l.i.start,
					  interval_l.i.end,
					  interval_l.i.strand,
					  false);		

	vector< UCSCElement<SV_BreakPoint*> > tmp_hits_r =
			r_bin.get(interval_r.i.chr,
					  interval_r.i.start,
					  interval_r.i.end,
					  interval_r.i.strand,
					  false);		

	if ( ( tmp_hits_l.size() == 0 ) ||
		 ( tmp_hits_r.size() == 0 ) ) {

		//insert(l_bin, r_bin);
		insert(l_bin, r_bin);

	} else { //Each side has at least one hit

		// Scan through the tmp hits to make sure the types match
		vector< UCSCElement<SV_BreakPoint*> >::iterator it;

		it = tmp_hits_l.begin();
		while(it != tmp_hits_l.end()) {
			if(it->value->type != type)
				it = tmp_hits_l.erase(it);
			else
				++it;
		}

		it = tmp_hits_r.begin();
		while(it != tmp_hits_r.end()) {
			if(it->value->type != type)
				it = tmp_hits_r.erase(it);
			else
				++it;
		}

		/*
		 * tmp_hits_left contains the pointers to the break points that have an
		 * interval overlapping the left end of the current pair,
		 * tmp_hits_right has the same for the right end.  To find the
		 * breakpoint(s) overlap both pairs we sort both vectors, then
		 * intersect the two vectors.  Since the vectors contain the address to
		 * the break point, the resulting intersection will contain a pointer
		 * to the object that needs to be modified. 
		 */

		sort(tmp_hits_l.begin(), tmp_hits_l.end(),
				UCSCElement<SV_BreakPoint*>::sort_ucscelement_by_value);

		sort(tmp_hits_r.begin(), tmp_hits_r.end(),
				UCSCElement<SV_BreakPoint*>::sort_ucscelement_by_value);

		vector< UCSCElement<SV_BreakPoint*> >::iterator t_it;

		vector< UCSCElement<SV_BreakPoint*> > 
			r(tmp_hits_l.size() + tmp_hits_r.size());

		it = set_intersection(tmp_hits_l.begin(), 
							  tmp_hits_l.end(),
							  tmp_hits_r.begin(),
							  tmp_hits_r.end(),
							  r.begin(),
							  UCSCElement<SV_BreakPoint*>::
									sort_ucscelement_by_value);

		if (it - r.begin() < 1) { //sides don't match, so insert

			insert(l_bin, r_bin);

		} else if (it - r.begin() == 1) { // both sides match same breakpoint
			// addr of the bp in both sets (and only item in the intersection)
			SV_BreakPoint *bp = r[0].value;

			vector< UCSCElement<SV_BreakPoint*> >::iterator l_it, r_it;

			/* 
			 * Merge the current break point with the one that was found to
			 * match (located at r[0] we need to remove r[0] from UCSCBin
			 * and re-add it with the new coordinates
			 */
			if ( bp->merge(this) ) {
				// find the left side match
				for (l_it = tmp_hits_l.begin(); l_it < tmp_hits_l.end(); ++l_it)
					if (l_it->value == bp)
						break;

				// find the right side match
				for (r_it = tmp_hits_r.begin(); r_it < tmp_hits_r.end(); ++r_it)
					if (r_it->value == bp)
						break;

				//int v = l_bin.remove(*l_it, false, false, true);
				//v = r_bin.remove(*r_it, false, false, true);
				l_bin.remove(*l_it, false, false, true);
				r_bin.remove(*r_it, false, false, true);

				delete this;

				bp->insert(l_bin, r_bin);
			} else {
			}
		}
	}
#endif
}
//}}}

//{{{ void SV_BreakPoint:: insert(UCSCBins<SV_BreakPoint*> &r_bin,
void
SV_BreakPoint::
insert(UCSCBins<SV_BreakPoint*> &r_bin)
{
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
