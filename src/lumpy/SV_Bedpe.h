/*****************************************************************************
 * SV_Bedpe.h
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

#ifndef __SV_BEDPE_H__
#define __SV_BEDPE_H__

#include "SV_Evidence.h"
#include "SV_BedpeReader.h"
#include "SV_BreakPoint.h"
#include "ucsc_bins.hpp"
#include "log_space.h"
#include "bedFilePE.h"

#include <iostream>

using namespace std;

class SV_Bedpe: public SV_Evidence
{
	friend ostream& operator<<(ostream& out, const SV_Bedpe& p);

	private:
		//void set_bp_interval_probability(struct breakpoint_interval *i);
		static void set_bp_interval_start_end(struct breakpoint_interval *i,
											  struct interval *target_interval,
											  struct interval *target_pair,
											  int back_distance,
											  int distro_size);
	public:
		SV_Bedpe(const BEDPE *bedpeEntry,
				 int weight,
				 int ev_id,
				 SV_BedpeReader *reader);
                void set_bp_interval_probability(struct breakpoint_interval *i);

		static void process_bedpe(const BEDPE *bedpeEntry,
                                          UCSCBins<SV_BreakPoint*> &r_bin,
                                          int weight,
                                          int ev_id,
                                          SV_BedpeReader *reader);


		static log_space* get_bp_interval_probability(char strand,
                                                              int distro_size,
                                                              double *distro);
		struct interval side_l;
		struct interval side_r;
		SV_BedpeReader *reader;

		bool is_interchromosomal();

		SV_BreakPoint* get_bp();

		void print_bedpe(int score);
		void print_evidence();
		string evidence_type();
};

#endif
