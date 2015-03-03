/*****************************************************************************
 * SV_BreakPoint.h
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

class SV_Evidence;
#ifndef __SV_BREAKPOINT_H__
#define __SV_BREAKPOINT_H__

#include "ucsc_bins.hpp"

#include "log_space.h"
#include "SV_Evidence.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

struct interval {
	CHR_POS start, end;
	CHR_POS start_clip, end_clip;
	char strand;
	string chr;
};

struct breakpoint_interval
{
	log_space *p;
	struct interval i;
};

class SV_BreakPoint
{
	friend ostream& operator<<(ostream& out, const SV_BreakPoint& b);

	private:
		static void trim_interval(log_space *interval_v,
			                  unsigned int size,
                                          int *trim_start,
                                          int *trim_end);

	public:
		static const int DELETION = 1;
		static const int DUPLICATION = 2;
		static const int INVERSION = 3;
		static const int TRANSLOCATION = 4;
		static double p_trim_threshold;
		static double p_merge_threshold;
		static bool sort_bp_by_interval_l_start(SV_BreakPoint *i,
                                                        SV_BreakPoint *j);
		static string ascii_interval_prob(
                        struct breakpoint_interval *i);
		static string ascii_prob(double *d, int size);
		int type;
		int weight;
		vector<SV_Evidence*> evidence;
		map<int, int> ev_ids;
		struct breakpoint_interval interval_l, interval_r;
		SV_BreakPoint(SV_Evidence *e);
		SV_BreakPoint();
		SV_BreakPoint(SV_BreakPoint *a, SV_BreakPoint *b);
		~SV_BreakPoint();
		void free_evidence();
                static void get_score( vector<SV_BreakPoint *> &bps,
                                       double *score_l,
                                       double *score_r);
                static void get_mixture(vector<SV_BreakPoint *> &bps,
                                        CHR_POS *start_l,
                                        CHR_POS *start_r,
                                        CHR_POS *end_l,
                                        CHR_POS *end_r,
                                        log_space **l,
                                        log_space **r);
                static int get_product(vector<SV_BreakPoint *> &bps,
                                        CHR_POS *start_l,
                                        CHR_POS *start_r,
                                        CHR_POS *end_l,
                                        CHR_POS *end_r,
                                        log_space **l,
                                        log_space **r);
                static bool does_intersect(struct breakpoint_interval *a,
                                           struct breakpoint_interval *b,
                                           bool check_strand);
                static void get_max_range(struct breakpoint_interval *b,
                                          CHR_POS *max_start,
                                          CHR_POS *max_end,
                                          log_space *max_value);
                bool merge(SV_BreakPoint *p);
                static bool test_interval_merge(
                        struct breakpoint_interval *curr_intr,
                        struct breakpoint_interval *new_intr,
                        CHR_POS *merged_start,
                        CHR_POS *merged_end,
                        log_space **merged_prob);
                void print_evidence(string pre);
                int trim_intervals();
                void init_interval_probabilities();
                void free_interval_probabilities();
                void print_bedpe(int id, int print_prob);
                void get_interval_probabilities(CHR_POS *start_l,
                                                CHR_POS *start_r,
                                                CHR_POS *end_l,
                                                CHR_POS *end_r,
                                                log_space **l,
                                                log_space **r);
                void print_interval_probabilities();
                void cluster(UCSCBins<SV_BreakPoint*> &r_bin);
                void insert(UCSCBins<SV_BreakPoint*> &l_bin,
                                        UCSCBins<SV_BreakPoint*> &r_bin);
                void insert(UCSCBins<SV_BreakPoint*> &r_bin);
                vector<int> get_evidence_ids();
                void do_it();
                void peak_distance(SV_BreakPoint *e,
                                   CHR_POS *dist,
                                   CHR_POS *overlap);
                static pair<CHR_POS,CHR_POS> min_pair(
                        vector< vector< pair<CHR_POS,CHR_POS> > > &m);
                void get_distances(vector< SV_BreakPoint*> &new_v);
                static void interval_product(struct breakpoint_interval *a_intr,
                                             struct breakpoint_interval *b_intr,
                                             CHR_POS *product_len,
                                             CHR_POS *product_start,
                                             CHR_POS *product_end,
                                             log_space **product_prob);
                void get_mixture(struct breakpoint_interval **l_mixture,
                                 struct breakpoint_interval **r_mixture);
		int get_max_sample_weight();
};
#endif
