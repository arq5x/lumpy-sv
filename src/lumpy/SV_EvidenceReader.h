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
#include <map>
#include "ucsc_bins.hpp"
#include "SV_BreakPoint.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

using namespace std;

class SV_EvidenceReader
{

public:
    static int counter;
    virtual ~SV_EvidenceReader();
    static map<int, string> sample_names;
    static map<int, string> ev_types;
    
    int ev_id;
    string sample_name;
    virtual string check_params();
    virtual bool add_param(char *param, char *val);

    virtual void initialize();
    virtual void set_statics();
    virtual void unset_statics();

    virtual void process_input_chr_pos(string chr,
                                       CHR_POS pos,
                                       UCSCBins<SV_BreakPoint*> &r_bin);
    virtual void process_input_chr_pos(string primary_chr,
                                       string secondary_chr,
                                       CHR_POS pos,
                                       UCSCBins<SV_BreakPoint*> &r_bin);
    virtual void process_input( BamAlignment &_bam,
                                RefVector &_refs,
                                BamWriter &inter_chrom_reads,
                                UCSCBins<SV_BreakPoint*> &r_bin);
    virtual void process_input( BamAlignment &_bam,
                                RefVector &_refs,
                                UCSCBins<SV_BreakPoint*> &r_bin);

    virtual void terminate();
    virtual string get_curr_chr();
    virtual CHR_POS get_curr_pos();

    virtual bool has_next();

    virtual string get_source_file_name();

    virtual string get_curr_primary_chr();
    virtual string get_curr_secondary_chr();

    virtual int32_t get_curr_primary_refid();
    virtual int32_t get_curr_secondary_refid();

    virtual CHR_POS get_curr_primary_pos();
    virtual CHR_POS get_curr_secondary_pos();

};
#endif
