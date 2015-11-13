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

#ifndef __SV_BAMREADER_H__
#define __SV_BAMREADER_H__

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

class SV_BamReader : public SV_EvidenceReader
{
    private:
        bool is_open, has_next_alignment;
	map<pair<string,string>, SV_EvidenceReader*> *bam_evidence_readers;
	BamMultiReader bam_reader;
	BamAlignment bam;
	string header;
	string bam_sort_order;
	SV_EvidenceReader *curr_reader;
	BamWriter inter_chrom_reads;
	string tmp_file_name;

    public:
        RefVector refs;
	SV_BamReader();
	~SV_BamReader();
	SV_BamReader(map<pair<string,string>,
                     SV_EvidenceReader*> *_bam_evidence_readers);

	bool add_param(char *param, char *val);
	string check_params();

	void initialize();
	void set_statics();

#if 0
	void process_input(UCSCBins<SV_BreakPoint*> &r_bin);
	void process_input(BamAlignment &_bam,
			   RefVector &_ref,
			   string header,
			   UCSCBins<SV_BreakPoint*> &r_bin);
        void process_input_chr(string chr,
                               UCSCBins<SV_BreakPoint*> &r_bin);
#endif
        void process_input_chr_pos(string chr,
                                   CHR_POS pos,
                                   UCSCBins<SV_BreakPoint*> &r_bin);

        void terminate();
        string get_curr_chr();
        CHR_POS get_curr_pos();
        bool has_next();
        string get_source_file_name();
        void set_inter_chrom_file_name(string file_name);
};
#endif
