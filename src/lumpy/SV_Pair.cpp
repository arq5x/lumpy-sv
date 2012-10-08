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
using namespace BamTools;

#include "SV_BreakPoint.h"
#include "SV_Pair.h"
#include "log_space.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <math.h>
#include <gsl_randist.h>
#include <gsl_cdf.h>
#include <gsl_statistics_double.h>

using namespace std;

//{{{ statics
//distro_sizedouble  SV_Pair:: cluster_Z = 0;
double  SV_Pair:: insert_Z = 0;
double  SV_Pair:: insert_stdev = 0;
double  SV_Pair:: insert_mean = 0;
int     SV_Pair:: min_non_overlap = 0;
double *SV_Pair:: histo = NULL;
double *SV_Pair:: distro = NULL;
int     SV_Pair:: distro_size = 0;
int     SV_Pair:: histo_size = 0;
int     SV_Pair:: histo_start = 0;
int     SV_Pair:: histo_end = 0;
int     SV_Pair:: back_distance = 0;
int     SV_Pair:: read_length = 0;
int     SV_Pair:: min_mapping_threshold = 0;
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
		int _id)
{
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

	if ( tmp_a.chr.compare(tmp_b.chr) > 0 ) {
		read_l = tmp_a;
		read_r = tmp_b;
	} else if ( tmp_a.chr.compare(tmp_b.chr) < 0 ) {
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
	id = _id;
}
//}}}

//{{{ void SV_Pair:: set_interval_probability()
void
SV_Pair::
set_bp_interval_probability(struct breakpoint_interval *i)
{
	int size = i->i.end - i->i.start + 1;
	log_space *tmp_p = (log_space *) malloc(size * sizeof(log_space));
	unsigned int j;
	for (j = 0; j < size; ++j) {
		if (i->i.strand == '+') 
			tmp_p[j] = get_ls(distro[j]);
		else
			tmp_p[(size - 1) - j] = get_ls(distro[j]);
	}

	i->p = tmp_p;
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
						  struct interval *target_pair)
{
#if 0
	cerr << target_interval->start << "\t" <<
			target_interval->end << "\t" <<
			target_pair->start << "\t" <<
			target_pair->end << endl;
#endif
			
	i->i.chr = target_interval->chr;
	i->i.strand = target_interval->strand;
	if ( i->i.strand == '+' ) {
		i->i.start = target_interval->end - back_distance;
		i->i.end = i->i.start + distro_size;
	} else {
		i->i.end = target_interval->start + back_distance;
		i->i.start = i->i.end - distro_size;
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

	set_bp_interval_start_end(&(new_bp->interval_l), &read_l, &read_r);
	set_bp_interval_start_end(&(new_bp->interval_r), &read_r, &read_l);

	set_bp_interval_probability(&(new_bp->interval_l));
	set_bp_interval_probability(&(new_bp->interval_r));

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

	if ( (read_r.end - read_l.start) >= insert_mean + (insert_Z*insert_stdev) )
		return true;

	if ( (read_r.end - read_l.start) <= insert_mean - (insert_Z*insert_stdev) )
		return true;

	return false;
}
//}}}

//{{{ bool SV_Pair:: is_sane()
bool
SV_Pair::
is_sane()
{
	if ( min_mapping_quality < min_mapping_threshold )
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


	int overlap = min(read_l.end, read_r.end) - max(read_l.start, read_r.start);
	//int non_overlap = min(read_len_a, read_len_b) - overlap;
	// how much does not overlap // overlap is at most read_len
	int non_overlap = min(read_len_a, read_len_b) - overlap;

	if ( (overlap > 0) && (abs(non_overlap) < min_non_overlap) )
		return false;
	else 
		return true;
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
#endif 

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
		"id:" << id << sep <<
		"weight:" << weight <<
		endl;
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

//{{{ void process_pair(const BamAlignment &curr,
void 
SV_Pair::
process_pair(const BamAlignment &curr,
			const RefVector refs,
			map<string, BamAlignment> &mapped_pairs,
			UCSCBins<SV_BreakPoint*> &l_bin,
			UCSCBins<SV_BreakPoint*> &r_bin,
			int weight,
			int id)
{
	if (mapped_pairs.find(curr.Name) == mapped_pairs.end())
		mapped_pairs[curr.Name] = curr;
	else {
		SV_Pair *new_pair = new SV_Pair(mapped_pairs[curr.Name],
										curr,
										refs,
										weight,
										id);
		if ( new_pair->is_sane() &&  new_pair->is_aberrant() ) {
			SV_BreakPoint *new_bp = new_pair->get_bp();
#ifdef TRACE
			cerr << "PE\t" << *new_bp << endl;
#endif
			new_bp->cluster(l_bin, r_bin);
		} else {
			free(new_pair);
		}

		mapped_pairs.erase(curr.Name);
	}
}
//}}}

#if 0
//{{{void SV_Pair:: cluster(UCSCBins<SV_BreakPoint*> &bins);
void
SV_Pair::
cluster(UCSCBins<SV_BreakPoint*> &l_bin,
		UCSCBins<SV_BreakPoint*> &r_bin)
{
	SV_BreakPoint *new_bp = get_bp();
	new_bp->cluster(l_bin, r_bin);

	// Find matching break points
	vector< UCSCElement<SV_BreakPoint*> > tmp_hits_a =
			bins.get(new_bp->interval_l.i.chr,
					 new_bp->interval_l.i.start,
					 new_bp->interval_l.i.end,
					 new_bp->interval_l.i.strand,
					 true);		

	vector< UCSCElement<SV_BreakPoint*> > tmp_hits_b =
			bins.get(new_bp->interval_r.i.chr,
					 new_bp->interval_r.i.start,
					 new_bp->interval_r.i.end,
					 new_bp->interval_r.i.strand,
					 true);		

	if ( ( tmp_hits_a.size() == 0 ) ||
		 ( tmp_hits_b.size() == 0 ) ) {
		cerr << "i1" << endl;
		insert(bins, new_bp);
	} else { //Each side has at least one hit
		/*
		 * tmp_hits_left contains the pointers to the break points that have an
		 * interval overlapping the left end of the current pair,
		 * tmp_hits_right has the same for the right end.  To find the
		 * breakpoint(s) overlap both pairs we sort both vectors, then
		 * intersect the two vectors.  Since the vectors contain the address to
		 * the break point, the resulting intersection will contain a pointer
		 * to the object that needs to be modified. 
		 */
		sort(tmp_hits_a.begin(), tmp_hits_a.end(),
				UCSCElement<SV_BreakPoint*>::sort_ucscelement_by_value);

		sort(tmp_hits_b.begin(), tmp_hits_b.end(),
				UCSCElement<SV_BreakPoint*>::sort_ucscelement_by_value);

		vector< UCSCElement<SV_BreakPoint*> > 
			r(tmp_hits_a.size() + tmp_hits_b.size());

		vector< UCSCElement<SV_BreakPoint*> >::iterator t_it;
		for (t_it = tmp_hits_a.begin(); t_it < tmp_hits_a.end(); ++t_it)
				cerr << "a:\t" << *(t_it->value) << endl;
		for (t_it = tmp_hits_b.begin(); t_it < tmp_hits_b.end(); ++t_it)
				cerr << "b:\t" << *(t_it->value) << endl;


		vector< UCSCElement<SV_BreakPoint*> >::iterator it;

		it = set_intersection(tmp_hits_a.begin(), 
							  tmp_hits_a.end(),
							  tmp_hits_b.begin(),
							  tmp_hits_b.end(),
							  r.begin(),
							  UCSCElement<SV_BreakPoint*>::
									sort_ucscelement_by_value);

		if (it - r.begin() < 1) { //sides don't match, so insert
			cerr << "i2" << endl;
			insert(bins, new_bp);
		} else if (it - r.begin() == 1) { // both sides match same breakpoint

			// addr of the bp in both sets (and only item in the intersection)
			SV_BreakPoint *bp = r[0].value;

			vector< UCSCElement<SV_BreakPoint*> >::iterator l_it, r_it;
			// find the left side match
			for (l_it = tmp_hits_a.begin();
					l_it < tmp_hits_a.end(); ++l_it)
				if (l_it->value == bp)
					break;

			// find the right side match
			for (r_it = tmp_hits_b.begin();
					r_it < tmp_hits_b.end(); ++r_it)
				if (r_it->value == bp)
					break;

			/* 
			 * Test to see if both sides hit the same interval 
			 * For this to be a valid hit, both ends need to hit the same
			 * breakpoint, but different ends of the breakpoint.  Each end has
			 * an id, test to see if the ids match
			 */
			if (l_it->id == r_it->id) { // both sides hit the same interval
				cerr << "i3" << endl;
				insert(bins, new_bp);
			}  else {
				cerr << "c" << endl;
				/* 
				 * Merge the current break point with the one that was found to
				 * match (located at r[0] we need to remove r[0] from UCSCBin
				 * and re-add it with the new coordinates
				 */
				if ( bp->merge(new_bp) ) {
					//cerr << bp << "\t" << *this << endl;
					//cerr << bp << "\t" << *bp << endl;
					//cerr << bp->evidence.size() << endl;
					free(new_bp);

					int v = bins.remove(*l_it, false, true, true);

					v = bins.remove(*r_it, false, true, true);

					insert(bins, bp);

				} else {
				}

			}
		}
	}
}
//}}}

//{{{ void SV_Pair:: insert(UCSCBins<SV_BreakPoint*> &bins, SV_BreakPoint
void
SV_Pair::
insert(UCSCBins<SV_BreakPoint*> &bins, SV_BreakPoint *new_bp)
{
	bins.add(new_bp->interval_l.i.chr,
			 new_bp->interval_l.i.start,
			 new_bp->interval_l.i.end,
			 new_bp->interval_l.i.strand,
			 new_bp);

	bins.add(new_bp->interval_r.i.chr,
			 new_bp->interval_r.i.start,
			 new_bp->interval_r.i.end,
			 new_bp->interval_r.i.strand,
			 new_bp);
}
//}}}
#endif

