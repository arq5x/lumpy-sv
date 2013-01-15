/*****************************************************************************
 * SV_EvidenceReader.h
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

#ifndef __SV_EVIDENCEREADER_H__
#define __SV_EVIDENCEREADER_H__

#include <string>
#include "ucsc_bins.hpp"
#include "SV_BreakPoint.h"
#include "api/BamReader.h"
using namespace BamTools;

using namespace std;

class SV_EvidenceReader
{

	public:
		static int counter;
		int sample_id;
		virtual string check_params();
		virtual bool add_param(char *param, char *val);

		virtual void initialize();
		virtual void set_statics();
		virtual void unset_statics();
		virtual void process_input( BamAlignment &_bam,
									RefVector &_ref,
									UCSCBins<SV_BreakPoint*> &r_bin);
		virtual void process_input( UCSCBins<SV_BreakPoint*> &r_bin);
		virtual void process_input_chr(string chr,
									   UCSCBins<SV_BreakPoint*> &r_bin);
		virtual void process_input_chr_pos(string chr,
										   CHR_POS pos,
										   UCSCBins<SV_BreakPoint*> &r_bin);
		virtual void terminate();
		virtual string get_curr_chr();
		virtual CHR_POS get_curr_pos();

		virtual bool has_next();

		virtual string get_source_file_name();
};
#endif
