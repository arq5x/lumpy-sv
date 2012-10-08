/*****************************************************************************
 * SV_BedpeReader.h
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

#ifndef __SV_BEDPEREADER_H__
#define __SV_BEDPEREADER_H__

#include <string>
using namespace std;

#include "ucsc_bins.hpp"
#include "SV_EvidenceReader.h"
#include "bedFilePE.h"

//{{{struct bedpe_parameters {
struct bedpe_parameters {
	string bedpe_file,
		   distro_file;
	unsigned int back_distance;
	int weight;
	int id;
};
//}}}

class SV_BedpeReader : public SV_EvidenceReader
{
	private:
		string bedpe_file,
			   distro_file;
		unsigned int back_distance;
		int weight;
		int id;
		bool is_open,
			 have_next_alignment;
	public:

	BedFilePE *bedpe;
	BEDPE bedpeEntry, nullBedpe;
	int lineNum;
	BedLineStatus bedpeStatus;

	SV_BedpeReader();
	bool add_param(char *param, char *val);
	string check_params();
	struct bedpe_parameters get_bedpe_parameters();

	void initialize();
	void set_statics();
	void process_input(UCSCBins<SV_BreakPoint*> &l_bin,
					   UCSCBins<SV_BreakPoint*> &r_bin);

	void process_input_chr(string chr,
						   UCSCBins<SV_BreakPoint*> &l_bin,
						   UCSCBins<SV_BreakPoint*> &r_bin);
	void terminate();
	string get_curr_chr();
	bool has_next();
};

#endif
