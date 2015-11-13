/*****************************************************************************
 * SV_SplitReadReader.h
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

#ifndef __SV_SPLITREADREADER_H__
#define __SV_SPLITREADREADER_H__

#include <string>
#include <map>
using namespace std;

#include "SV_BreakPoint.h"
#include "SV_EvidenceReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "ucsc_bins.hpp"

using namespace BamTools;

//{{{ struct split_read_parameters {
struct split_read_parameters {
	string bam_file,
		sample_name;
	unsigned int min_non_overlap,
				 back_distance,
				 min_mapping_threshold,
                 min_clip;
	int weight;
	vector<string> read_group;
};
//}}}

class SV_SplitReadReader : public SV_EvidenceReader
{
    public:
        string bam_file;
        unsigned int min_non_overlap,
                     back_distance,
                     min_mapping_threshold,
                 min_clip;
        int weight;
        vector<string> read_group;
        bool is_open,
        have_next_alignment;

        BamAlignment bam;
        BamReader reader;
        map<string, BamAlignment> mapped_splits;
        string header;
        RefVector refs;
        bool inited;

        ~SV_SplitReadReader();
        SV_SplitReadReader();
        bool add_param(char *param, char *val);
        string check_params();
        void initialize();
        void set_statics();
        void unset_statics();
        //void process_input( UCSCBins<SV_BreakPoint*> &r_bin);

        void process_input( BamAlignment &_bam,
        RefVector &_ref,
        UCSCBins<SV_BreakPoint*> &r_bin);

        void process_input(BamAlignment &_bam,
                           RefVector &_ref,
                           BamWriter &inter_chrom_reads,
                           UCSCBins<SV_BreakPoint*> &r_bin);

        //void process_input_chr(string chr,
        //UCSCBins<SV_BreakPoint*> &r_bin);
        //void process_input_chr_pos(string chr,
        //CHR_POS pos,
        //UCSCBins<SV_BreakPoint*> &r_bin);
        void terminate();
        string get_curr_chr();
        CHR_POS get_curr_pos();
        bool has_next();
        string get_source_file_name();
};

#endif
