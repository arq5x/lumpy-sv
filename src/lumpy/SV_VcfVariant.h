/*****************************************************************************
 * SV_Vcf.h
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

#include "SV_Vcf.h"
#include "SV_BreakPoint.h"

#include <string>
#include <vector>
using namespace std;

class SV_VcfVariant
{
 public:
    
    SV_VcfVariant(SV_Vcf *vcf);
    SV_VcfVariant(SV_Vcf *vcf,
		  SV_BreakPoint *bp,
		  int bp_id,
		  int print_prob);
    ~SV_VcfVariant();

    SV_Vcf *vcf;

    // variant fields
    string chrom;
    uint64_t pos;
    string id;
    string ref;
    string alt;
    string qual; // string to allow for '.'
    string filter;
    string info;
    vector<string> format;
    
    void add_info(string id);
    void add_info(string id,
		  string value);

    void add_format_field(string fmt);

    void print_var();
    
};

#endif
