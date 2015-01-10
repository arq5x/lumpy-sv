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
add_info(string id)
{
}

void
SV_Vcf::
add_info(string id,
	 string value)
{
}


     

void
SV_Vcf::
print_header()
{
    string sep = "\t";

    cout << "##fileformat=VCFv4.2" << endl <<
	// "##INFO=<ID=TOOL,Number=1,Type=String,Description=\"Tool used to generate variant call\">" << endl <<
	"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl <<
	"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl <<
	"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl <<
	// "##INFO=<ID=STR,Number=.,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format\">" << endl <<
	"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << endl <<
	"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << endl <<
	"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << endl <<
	// "##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">" << endl <<
	// "##INFO=<ID=PARID,Number=1,Type=String,Description=\"ID of partner breakend\">" << endl <<
	"##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">" << endl <<
	"##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">" << endl <<
	// "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">" << endl <<
	// "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">" << endl <<
	// "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">" << endl <<
	// "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">" << endl <<
	// "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency, for each ALT allele, in the same order as listed\">" << endl <<
	// "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">" << endl <<
	// "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl <<
	"##INFO=<ID=SUP,Number=.,Type=Integer,Description=\"Number of pieces of evidence supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=PESUP,Number=.,Type=Integer,Description=\"Number of paired-end reads supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=SRSUP,Number=.,Type=Integer,Description=\"Number of split reads supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=EVTYPE,Number=.,Type=String,Description=\"Type of LUMPY evidence contributing to the variant call\">" << endl <<
	"##INFO=<ID=PRIN,Number=0,Type=Flag,Description=\"Indicates variant as the principal variant in a BEDPE pair\">" << endl <<
	"##ALT=<ID=DEL,Description=\"Deletion\">" << endl <<
	"##ALT=<ID=DUP,Description=\"Duplication\">" << endl <<
	"##ALT=<ID=INV,Description=\"Inversion\">" << endl <<
	"##ALT=<ID=DUP:TANDEM,Description=\"Tandem duplication\">" << endl <<
	"##ALT=<ID=INS,Description=\"Insertion of novel sequence\">" << endl <<
	"##ALT=<ID=CNV,Description=\"Copy number variable region\">" << endl <<
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl <<
	"##FORMAT=<ID=SUP,Number=1,Type=Integer,Description=\"Number of pieces of evidence supporting the variant\">" << endl <<
	"##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Number of paired-end reads supporting the variant\">" << endl <<
	"##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of split reads supporting the variant\">" << endl <<
	"##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality\">" << endl;
    // "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">" << endl <<
    // "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">" << endl <<
    // "##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">" << endl <<
    // "##FORMAT=<ID=CNL,Number=.,Type=Float,Description=\"Copy number genotype likelihood form imprecise events\">" << endl <<
    // "##FORMAT=<ID=NQ,Number=1,Type=Integer,Description=\"Phred style probability score that the variant is novel\">" << endl <<
    // "##FORMAT=<ID=HAP,Number=1,Type=Integer,Description=\"Unique haplotype identifier\">" << endl <<
    // "##FORMAT=<ID=AHAP,Number=1,Type=Integer,Description=\"Unique identifier of ancestral haplotype\">" << endl <<
    
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
