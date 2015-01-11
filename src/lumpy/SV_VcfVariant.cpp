/*****************************************************************************
 * SV_VcfVariant.cpp
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

#include "SV_VcfVariant.h"
#include "SV_BreakPoint.h"
#include "SV_Evidence.h"

#include <vector>
#include <string>
#include <sstream>
using namespace std;

SV_VcfVariant::
SV_VcfVariant(SV_Vcf *_vcf)
{
    vcf = _vcf;
    chrom = ".";
    pos = 0;
    id = ".";
    ref = "N";
    alt = ".";
    qual = ".";
    filter = ".";
    info = "";
}

SV_VcfVariant::
SV_VcfVariant(SV_Vcf *_vcf,
	      SV_BreakPoint *bp,
	      int bp_id,
	      int print_prob)
{
    // set as object variable
    vcf = _vcf;
    
    map<string,int> uniq_strands;
    vector<SV_Evidence*>::iterator it;
    vector<SV_BreakPoint *> bps;

    for (it = bp->evidence.begin(); it < bp->evidence.end(); ++it) {
	SV_Evidence *e = *it;
	SV_BreakPoint *tmp_bp = e->get_bp();
	tmp_bp->init_interval_probabilities();

	bps.push_back(tmp_bp);

	stringstream strands;
	strands << tmp_bp->interval_l.i.strand << tmp_bp->interval_r.i.strand;

	if (uniq_strands.find(strands.str()) == uniq_strands.end())
	    uniq_strands[strands.str()] = 1;
	else
	    uniq_strands[strands.str()] = uniq_strands[strands.str()] + 1;
    }

    double score_l, score_r;
    bp->get_score(bps, &score_l, &score_r);

    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
	SV_BreakPoint *tmp_bp = *bp_it;
	delete tmp_bp;
    }

    // Get the most likely posistions
    CHR_POS start_l, start_r, end_l, end_r;
    log_space *l,*r;
    bp->get_interval_probabilities(&start_l,
			       &start_r,
			       &end_l,
			       &end_r,
			       &l,
			       &r);

    // l and r contain the full, untrimmed distribution.  All of the
    // calculations below must be shifted to the right by an offset so they are
    // considering the correction area of the pdf
    CHR_POS l_trim_offset = bp->interval_l.i.start - start_l;
    CHR_POS r_trim_offset = bp->interval_r.i.start - start_r;

    // find relative position of max value in left
    log_space max = -INFINITY;
    unsigned int i, l_max_i = 0, r_max_i = 0;
    for (i = 0; i < (bp->interval_l.i.end - bp->interval_l.i.start + 1); ++i) {
	if (l[l_trim_offset + i] > max) {
	    max = l[l_trim_offset + i];
	    l_max_i = l_trim_offset + i;
	}
    }

    // need to use start_l here, not interval_l.i.start, since l[0] is at
    // position start_l, and l_max_i is w.r.t. l[0]
    CHR_POS abs_max_l = start_l + l_max_i;

    // find relative position of max value in right
    max = -INFINITY;
    for (i = 0; i < (bp->interval_r.i.end - bp->interval_r.i.start + 1); ++i) {
	if (r[r_trim_offset + i] > max) {
	    max = r[r_trim_offset + i];
	    r_max_i = r_trim_offset + i;
	}
    }

    // need to use start_r here, not interval_r.i.start, since r[0] is at
    // position start_r, and r_max_i is w.r.t. r[0]
    CHR_POS abs_max_r = start_r + r_max_i;

    // stream to convert ints to strings
    ostringstream convert;

    // set first 7 VCF columns
    chrom = bp->interval_l.i.chr;
    pos = abs_max_l;
    convert.str(string());
    convert << bp_id;
    id = convert.str(); // Note: removed the id: -1 control statement
    ref = "N";
    if (bp->interval_l.i.chr.compare(bp->interval_r.i.chr) != 0) {
	alt = "INTERCHROM";
	add_info("SVTYPE", "BND");
    }
    else if (bp->type == bp->DELETION ) {
	alt = "<DEL>";
	add_info("SVTYPE", "DEL");
    }
    else if (bp->type == bp->DUPLICATION) {
	alt = "<DUP>";
	add_info("SVTYPE", "DUP");
    }
    else if (bp->type == bp->INVERSION) {
	alt = "<INV>";
	add_info("SVTYPE", "INV");
    }
    else
	alt = ".";
    qual = ".";
    filter = ".";

    // INFO: SVLEN
    int64_t svlen = abs_max_l - abs_max_r;
    convert.str(string());
    convert << svlen;
    add_info("SVLEN", convert.str());

    // INFO: END
    convert.str(string());
    convert << abs_max_r;
    add_info("END", convert.str());

    // INFO: LP, RP
    if (print_prob > 0) {
	string lp, rp;
	for (i = 0;
	     i < (bp->interval_l.i.end - bp->interval_l.i.start + 1);
	     ++i) {
	    if (i != 0)
		lp.append(",");
	    convert.str(string());
	    convert << get_p(l[l_trim_offset+i]);
	    lp.append(convert.str());
	}
	for (i = 0;
	     i < (bp->interval_r.i.end - bp->interval_r.i.start + 1);
	     ++i) {
	    if (i != 0)
		rp.append(",");
	    convert.str(string());
	    convert << get_p(r[r_trim_offset+i]);
	    rp.append(convert.str());
	}
	add_info("LP", lp);
	add_info("RP", rp);
    }

    // FORMAT
    vector<string>::iterator samp_itr;
    for (samp_itr = vcf->samples.begin();
    	 samp_itr != vcf->samples.end();
    	 ++samp_itr) {
    	add_format(*samp_itr, "GT", "./.");
	add_format(*samp_itr, "SUP", "20");
    }
}

void
SV_VcfVariant::
add_format(string sample,
	   string format,
	   string value)
{
    // add to VCF's active format vector
    this->vcf->add_format(format);

    var_samples[sample][format] = value;
    
    // cout << sample << format << value << endl;
}

void
SV_VcfVariant::
add_info(string id)
{
    if (info != "")
	info.append(";");
    info.append("id");
}

void
SV_VcfVariant::
add_info(string id,
	 string val)
{
    if (info != "")
	info.append(";");
    info.append(id);
    info.append("=");
    info.append(val);
}

string
SV_VcfVariant::
get_format_string()
{
    ostringstream fmt_ss;
    vector<string>::iterator fmt_it;
    for (fmt_it = this->vcf->active_formats.begin();
	 fmt_it != this->vcf->active_formats.end();
	 ++fmt_it) {
	if (fmt_ss.str() != "")
	    fmt_ss << ":";
	fmt_ss  << *fmt_it;
    }

    return fmt_ss.str();
}

string
SV_VcfVariant::
get_sample_string(string sample)
{
    ostringstream samp_ss;
    vector<string>::iterator fmt_it;
    for (fmt_it = this->vcf->active_formats.begin();
	 fmt_it != this->vcf->active_formats.end();
	 ++fmt_it) {
	if (samp_ss.str() != "")
	    samp_ss << ":";
	samp_ss << var_samples[sample][*fmt_it];
    }

    return samp_ss.str();
}

void
SV_VcfVariant::
print_var()
{
    string sep = "\t";
    string format = get_format_string();
    
    cout << chrom << sep <<
	pos << sep <<
	id << sep <<
	ref << sep <<
	alt << sep <<
	qual << sep <<
	filter << sep <<
	info << sep <<
	format;
    
    vector<string>::iterator samp_it;
    for (samp_it = this->vcf->samples.begin();
	 samp_it != this->vcf->samples.end();
	 ++samp_it)
	cout << sep <<
	    get_sample_string(*samp_it);
    
    cout << endl;
}

