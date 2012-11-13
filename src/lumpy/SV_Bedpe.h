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
											  struct interval *target_pair);
	public:
		SV_Bedpe(const BEDPE *bedpeEntry,
				 int weight,
				 int id,
				 int sample_id);

		static void process_bedpe(const BEDPE *bedpeEntry,
								  UCSCBins<SV_BreakPoint*> &l_bin,
								  UCSCBins<SV_BreakPoint*> &r_bin,
								  int weight,
								  int id,
								  int sample_id);


		static log_space* get_bp_interval_probability(char strand);

		static double *distro;
		static int distro_size;
		static int distro_start;
		static int distro_end;
		static int back_distance;

		struct interval side_l;
		struct interval side_r;

		SV_BreakPoint* get_bp();

		void print_bedpe(int score);
		void print_evidence();
};

#endif
