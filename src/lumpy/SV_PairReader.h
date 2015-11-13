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

#ifndef __SV_PAIRREADER_H__
#define __SV_PAIRREADER_H__

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

//{{{struct pair_end_parameters {
struct pair_end_parameters {
    string bam_file,
           sample_name,
           histo_file;
    double mean, stdev;
    unsigned int read_length,
                 min_non_overlap,
                 discordant_z,
                 back_distance,
                 min_mapping_threshold;
    int weight;
    vector<string> read_group;
};
//}}}

class SV_PairReader : public SV_EvidenceReader
{
    public:
        string bam_file,
               histo_file;
        double mean, stdev;
        unsigned int read_length,
                     min_non_overlap,
                     discordant_z,
                     back_distance,
                     min_mapping_threshold;
        int weight;
        vector<string> read_group;
        bool is_open,
             have_next_alignment;
        double *histo;
        double *distro;
        unsigned int histo_start,histo_end;
        int distro_size;

        BamAlignment bam;
        BamReader reader;
        map<string, BamAlignment> mapped_pairs;
        string header;
        RefVector refs;
        bool inited;

        ~SV_PairReader();
        SV_PairReader();
        SV_PairReader(struct pair_end_parameters pair_end_param);
        bool add_param(char *param, char *val);
        string check_params();
        struct pair_end_parameters get_pair_end_parameters();
        void initialize();
        void set_statics();
        void unset_statics();
        //void process_input( UCSCBins<SV_BreakPoint*> &r_bin);
            
        void process_input(BamAlignment &_bam,
                           RefVector &_ref,
                           BamWriter &inter_chrom_reads,
                           UCSCBins<SV_BreakPoint*> &r_bin);

        void process_input(BamAlignment &_bam,
                           RefVector &_ref,
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
