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
    ~SV_VcfVariant();

    // variant fields
    string chrom;
    uint64_t pos;
    string id;
    string ref;
    string alt;
    string qual; // string to allow for '.'
    string filter;
    string info;

    void add_sample(string sample_name);
    void add_info(string id);
    void add_info(string id,
		  string value);
    void add_format(string sample,
		    string format,
		    string value);

    string get_format_string();
    string get_sample_string(string sample);

    void print_header();
    void print_var();

    // should be private
    vector<string> samples;
    vector<string> active_formats;
    map< string,map<string,string> > var_samples;

 /* private: */
};


#endif
