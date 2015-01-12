/*****************************************************************************
 * SV_VcfVariant.h
 * (c) 2014 - Ryan M. Layer
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

#ifndef __SV_VCFVARIANT_H__
#define __SV_VCFVARIANT_H__

#include "SV_BreakPoint.h"

#include <string>
#include <vector>
#include <map>
using namespace std;

class SV_VcfVariant
{
 public:
    SV_VcfVariant();
    SV_VcfVariant(SV_BreakPoint *bp,
		  int bp_id,
		  int print_prob);
    /* ~SV_VcfVariant(); */

    // variant fields
    string chrom;
    uint64_t pos;
    string id;
    string ref;
    string alt;
    string qual; // string to allow for '.'
    string filter;
    string info;

    map<string,int> get_strands(SV_BreakPoint *bp);
    vector<string> get_alt(string l_chrom,
			   CHR_POS l_pos,
			   string r_chrom,
			   CHR_POS r_pos,
			   string ref,
			   string strands);
    
    void add_sample(string sample_name);
    void set_info(string id);
    void set_info(string id,
		  string value);

    void set_sample_field(string sample,
			  string format,
			  string value);
    string get_sample_field(string sample,
			    string format);
    
    string get_format_string();
    string get_sample_string(string sample);

    void print_header();
    void print_var();

 private:
    static vector<string> samples;
    static vector<string> active_formats;
    map< string,map<string,string> > sample_vals;
};


#endif
