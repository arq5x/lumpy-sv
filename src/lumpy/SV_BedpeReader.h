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
	    sample_name;
        /*
	string bedpe_file,
		   distro_file;
        */
	unsigned int back_distance;
	int weight;
};
//}}}

class SV_BedpeReader : public SV_EvidenceReader
{
	public:
		string bedpe_file;
                /* 
		string bedpe_file,
			   distro_file;
                */
		unsigned int back_distance;
		int weight;
		bool is_open,
			 have_next_alignment;

		double *histo;
		int histo_start,histo_end;
		log_space *distro;
		int distro_size;
		int distro_start;
		int distro_end;


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
	void unset_statics();
	void process_input( UCSCBins<SV_BreakPoint*> &r_bin);

	void process_input_chr(string chr,
						   UCSCBins<SV_BreakPoint*> &r_bin);
	void process_input_chr_pos(string chr,
							   CHR_POS pos,
							   UCSCBins<SV_BreakPoint*> &r_bin);
	void terminate();
	string get_curr_chr();
	CHR_POS get_curr_pos();
	bool has_next();

};

#endif
