/*****************************************************************************
 * SV_PairReader.h
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

#ifndef __SV_INTERCHROMBAMREADER_H__
#define __SV_INTERCHROMBAMREADER_H__

#include <string>
#include <map>
using namespace std;

#include "SV_BreakPoint.h"
#include "SV_EvidenceReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamMultiReader.h"
#include "api/BamAux.h"

#include "ucsc_bins.hpp"

using namespace BamTools;

class SV_InterChromBamReader : public SV_EvidenceReader
{
	private:
		bool is_open, has_next_alignment;
		map<pair<string,string>, SV_EvidenceReader*> *bam_evidence_readers;
		BamReader bam_reader;
		RefVector refs;
		BamAlignment bam;
		string header;
		string inter_chrom_bam_file;

	public:
		string get_curr_primary_chr();
		string get_curr_secondary_chr();

		int32_t get_curr_primary_refid();
		int32_t get_curr_secondary_refid();

		CHR_POS get_curr_primary_pos();
		CHR_POS get_curr_secondary_pos();

		SV_InterChromBamReader();
		~SV_InterChromBamReader();
		SV_InterChromBamReader(
				string _inter_chrom_bam_file,
				map<pair<string,string>, SV_EvidenceReader*> *_bam_evidence_readers);

		void process_input_chr_pos(string primary_chr,
								   string secondary_chr,
								   CHR_POS pos,
								   UCSCBins<SV_BreakPoint*> &r_bin);
		bool add_param(char *param, char *val);
		string check_params();

		void initialize();
		void set_statics();

#if 0
		void process_input( UCSCBins<SV_BreakPoint*> &r_bin);
		void process_input( BamAlignment &_bam,
							RefVector &_ref,
							string header,
							UCSCBins<SV_BreakPoint*> &r_bin);

		void process_input_chr(string chr,
							   UCSCBins<SV_BreakPoint*> &r_bin);
		void process_input_chr_pos(string chr,
								   CHR_POS pos,
								   UCSCBins<SV_BreakPoint*> &r_bin);
#endif

		void terminate();
		string get_curr_chr();
		//CHR_POS get_curr_pos();
		bool has_next();
		string get_source_file_name();
};

#endif
