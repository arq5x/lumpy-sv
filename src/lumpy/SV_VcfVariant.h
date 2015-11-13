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
    ~SV_VcfVariant();
    SV_VcfVariant();
    SV_VcfVariant(SV_BreakPoint *bp,
		  int bp_id,
		  int print_prob);
    /* ~SV_VcfVariant(); */

    bool is_multiline;
    
    // variant fields (array for multi-line vars)
    string chrom[2];
    CHR_POS pos[2];
    string id[2];
    string ref[2];
    string alt[2];
    string qual[2]; // string to allow for '.'
    string filter[2];
    string info[2];

    map<string,int> get_strands(SV_BreakPoint *bp);
    vector<string> get_alt(string l_chrom,
			   CHR_POS l_pos,
			   string l_ref,
			   string r_chrom,
			   CHR_POS r_pos,
			   string r_ref,
			   string strands);
    
    void add_sample(string sample_name);
    void set_info(int line,
		  string id);
    void set_info(int line,
		  string id,
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

    template <typename Type>
      string to_string(Type t);
    
 private:
    static const int LINE1 = 0;
    static const int LINE2 = 1;
    static vector<string> samples;
    static vector<string> active_formats;
    map< string,map<string,string> > sample_vals;
};


#endif
