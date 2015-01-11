/*****************************************************************************
 * SV_Vcf.cpp
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

#include "SV_Vcf.h"

#include <string>
#include <vector>
#include <iostream>
using namespace std;

// SV_Vcf::
// SV_Vcf() {}

void
SV_Vcf::
add_sample(string sample_name)
{
    if (find(samples.begin(),
	     samples.end(),
	     sample_name)
	== samples.end())
	samples.push_back(sample_name);
}

void
SV_Vcf::
add_format(string format)
{
    if (find(active_formats.begin(),
	     active_formats.end(),
	     format)
	== active_formats.end())
	active_formats.push_back(format);
}

void
SV_Vcf::
print_header()
{
    string sep = "\t";

    cout << "##fileformat=VCFv4.2" << endl <<
	"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl <<
	"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl <<
	"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl <<
	"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << endl <<
	"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << endl <<
	"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << endl <<
	"##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">" << endl <<
	"##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">" << endl <<
	"##INFO=<ID=SP,Number=.,Type=Integer,Description=\"Number of pieces of evidence supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=PE,Number=.,Type=Integer,Description=\"Number of paired-end reads supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=SR,Number=.,Type=Integer,Description=\"Number of split reads supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=EV,Number=.,Type=String,Description=\"Type of LUMPY evidence contributing to the variant call\">" << endl <<
	"##INFO=<ID=LP,Number=.,Type=String,Description=\"LUMPY probability curve of the left breakend\">" << endl <<
	"##INFO=<ID=RP,Number=.,Type=String,Description=\"LUMPY probability curve of the right breakend\">" << endl <<
	// "##INFO=<ID=PRIN,Number=0,Type=Flag,Description=\"Indicates variant as the principal variant in a BEDPE pair\">" << endl <<
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl <<
	"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Number of pieces of evidence supporting the variant\">" << endl <<
	"##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Number of paired-end reads supporting the variant\">" << endl <<
	"##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of split reads supporting the variant\">" << endl;
    
    cout <<
	"#CHROM" << sep <<
	"POS" << sep <<
	"ID" << sep <<
	"REF" << sep <<
	"ALT" << sep <<
	"QUAL" << sep <<
	"FILTER" << sep <<
	"INFO" << sep <<
	"FORMAT";

    for (vector<string>::iterator samp = samples.begin();
	 samp != samples.end();
	 ++samp) {
	cout << sep <<
	    *samp;
    }
    cout << endl;
}
