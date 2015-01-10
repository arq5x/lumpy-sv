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

#ifndef __SV_VCF_H__
#define __SV_VCF_H__

#include <string>
#include <vector>
using namespace std;

class SV_Vcf
{
 public:
    vector<string> samples;
    vector<string> active_formats;

    void add_sample(string sample_name);
    void print_header();
    void print_variant();

};

/* class SV_VcfHeader */
/* { */
/*  public: */
/*     VcfHeader(); */
/*     ~VcfHeader() {} */

/* } */

/* class SV_VcfVariant */
/* { */

/* } */

#endif
