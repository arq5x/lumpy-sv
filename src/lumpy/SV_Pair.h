/*****************************************************************************
 * SV_Pair.h
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

#ifndef __SV_PAIR_H__
#define __SV_PAIR_H__

//#include "BamAncillary.h"
using namespace BamTools;

#include "SV_Evidence.h"
#include "SV_PairReader.h"
#include "SV_BreakPoint.h"
#include "SV_PairReader.h"
#include "ucsc_bins.hpp"
#include "log_space.h"

#include <iostream>

using namespace std;

class SV_Pair: public SV_Evidence
{
    friend ostream& operator<<(ostream& out, const SV_Pair& p);

    private:
        //void set_bp_interval_probability(struct breakpoint_interval *i);
        static void set_bp_interval_start_end(struct breakpoint_interval *i,
                                              struct interval *target_interval,
                                              struct interval *target_pair,
                                              unsigned int back_distance,
                                              unsigned int distro_size);

	public:
		//static double insert_mean;
		//static double insert_stdev;
		static double insert_Z;
		//static int min_non_overlap;
		//static log_space *distro;
		//static double *histo;
		//static int distro_size;
		//static int histo_size;
		//static int histo_start;
		//static int histo_end;
		//static int back_distance;
		//static int read_length;
		//static int min_mapping_threshold;

		unsigned int min_mapping_quality;
		struct interval read_l;
		struct interval read_r;
		bool read_l_is_split, read_r_is_split;
		SV_PairReader *reader;


		SV_Pair(const BamAlignment &bam_a,
				const BamAlignment &bam_b,
				const RefVector &refs,
				int weight,
				int ev_id,
				SV_PairReader *reader);

		static void process_pair(const BamAlignment &curr,
								const RefVector refs,
								map<string, BamAlignment> &mapped_pairs,
								UCSCBins<SV_BreakPoint*> &r_bin,
								int weight,
								int ev_id,
								SV_PairReader *reader);

		static void process_intra_chrom_pair(
								 const BamAlignment &curr,
								 const RefVector refs,
								 BamWriter &inter_chrom_reads,
								 map<string, BamAlignment> &mapped_pairs,
								 UCSCBins<SV_BreakPoint*> &r_bin,
								 int weight,
								 int ev_id,
								 SV_PairReader *reader);

		static log_space* get_bp_interval_probability(char strand,
													  int distro_size,
													  double *distro);
		SV_BreakPoint* get_bp();

		bool is_aberrant();
		bool is_sane();
		bool is_interchromosomal();
		void cluster( UCSCBins<SV_BreakPoint*> &r_bin);

		static void set_distro_from_histo ();
		static int set_distro_from_histo(int back_distance,
										 int histo_start,
										 int histo_end,
										 double *histo,
										 double **distro);

		void print_evidence();

		void print_bedpe(int score);
		string evidence_type();
};

#endif
