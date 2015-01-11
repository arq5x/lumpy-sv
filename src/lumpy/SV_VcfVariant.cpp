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

// initialize static members
vector<string> SV_VcfVariant::samples;
vector<string> SV_VcfVariant::active_formats;

SV_VcfVariant::
SV_VcfVariant()
{
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
SV_VcfVariant(SV_BreakPoint *bp,
	      int bp_id,
	      int print_prob)
{
    map<string,int> uniq_strands;
    vector<SV_Evidence*>::iterator it;
    vector<SV_BreakPoint *> bps;
    vector<string>::iterator samp_it;

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

    // Get the most likely positions
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

    // set fixed VCF columns
    chrom = bp->interval_l.i.chr;
    pos = abs_max_l;
    id = to_string(bp_id); // Note: removed the id: -1 control statement
    ref = "N";
    if (bp->interval_l.i.chr.compare(bp->interval_r.i.chr) != 0) {
	alt = "INTERCHROM";
	set_info("SVTYPE", "BND");
    }
    else if (bp->type == bp->DELETION ) {
	alt = "<DEL>";
	set_info("SVTYPE", "DEL");
    }
    else if (bp->type == bp->DUPLICATION) {
	alt = "<DUP>";
	set_info("SVTYPE", "DUP");
    }
    else if (bp->type == bp->INVERSION) {
	alt = "<INV>";
	set_info("SVTYPE", "INV");
    }
    else
	alt = ".";
    qual = ".";
    filter = ".";

    // INFO: SVLEN
    int64_t svlen = abs_max_l - abs_max_r;
    set_info("SVLEN", to_string(svlen));

    // INFO: END
    set_info("END", to_string(abs_max_r));

    // FORMAT: initialize appropriate fields
    for (samp_it = samples.begin();
    	 samp_it != samples.end();
    	 ++samp_it) {
    	set_sample_field(*samp_it, "GT", "./.");
	set_sample_field(*samp_it, "SU", "0");
	string ev_arr[] = {"PE", "SR", "BD"};
	for (int i = 0; i < 3; ++i) {
	    if (find(active_formats.begin(),
		     active_formats.end(),
		     ev_arr[i])
		!= active_formats.end())
		set_sample_field(*samp_it, ev_arr[i], "0");
	}
    }

    // FORMAT: update variant support
    map<string,int> var_supp;
    map<int, int>::iterator ids_it;
    vector<int> _ids;
    for (ids_it = bp->ids.begin();
	 ids_it != bp->ids.end();
	 ++ids_it)
    	_ids.push_back(ids_it->first);
    vector<int>::iterator _ids_it;
    for (_ids_it = _ids.begin();
	 _ids_it != _ids.end();
	 ++_ids_it) {
	string samp = SV_EvidenceReader::sample_names[*_ids_it];
	string ev_type = SV_EvidenceReader::ev_types[*_ids_it];

	int samp_supp = atoi(get_sample_field(samp, ev_type).c_str());
	int ev = atoi(get_sample_field(samp, ev_type).c_str());
	int new_ev = bp->ids[*_ids_it];
	
	ev += new_ev;
	samp_supp += new_ev;

	var_supp["SU"] += new_ev;
	var_supp[ev_type] += new_ev;
	
	set_sample_field(samp, ev_type, to_string(ev));
	set_sample_field(samp, "SU", to_string(samp_supp));
    }
    
    // INFO: support
    map<string,int>::iterator ev_it;
    for (ev_it = var_supp.begin();
	 ev_it != var_supp.end();
	 ++ev_it) {
	set_info(ev_it->first, to_string(ev_it->second));
    }
    
    // += bp->ids[*_ids_it];
    // cout << "id: " << *_ids_it << endl;
    // cout << *_ids_it << ":" << bp->ids[*_ids_it] << endl;
    // set_sample_field(*samp_it, "SU", to_string(samp_support));
    // }

    // }
    

    // cout << "\t";

    // cout << "STRANDS:";
    // map<string,int>:: iterator s_it;
    // for ( s_it = uniq_strands.begin(); s_it != uniq_strands.end(); ++s_it) {
    // 	if (s_it != uniq_strands.begin())
    // 	    cout <<  ";";
    // 	cout << s_it->first << "," << s_it->second;
    // }

    // cout << "\t";

    // INFO: LP, RP
    if (print_prob > 0) {
	string lp, rp;
	for (i = 0;
	     i < (bp->interval_l.i.end - bp->interval_l.i.start + 1);
	     ++i) {
	    if (i != 0)
		lp.append(",");
	    lp.append(to_string(get_p(l[l_trim_offset+i])));
	}
	for (i = 0;
	     i < (bp->interval_r.i.end - bp->interval_r.i.start + 1);
	     ++i) {
	    if (i != 0)
		rp.append(",");
	    rp.append(to_string(get_p(r[r_trim_offset+i])));
	}
	set_info("LP", lp);
	set_info("RP", rp);
    }

    // INFO: IMPRECISE
    set_info("IMPRECISE"); // add conditional here

    // INFO: Get the area that includes the max and 95% of the probabitliy
    log_space p_95 = get_ls(0.95);
    log_space total = l[l_max_i];
    CHR_POS l_l_i = l_max_i,
	l_r_i = l_max_i,
	t_last = bp->interval_l.i.end - bp->interval_l.i.start;

    while ( ((l_l_i > 0) || (l_r_i < t_last)) && (total < p_95) ){
	if ( l_l_i == 0 ) {
	    total = ls_add(total, l[l_r_i+1]);;
	    ++l_r_i;
	} else if ( l_r_i == t_last ) {
	    total = ls_add(total, l[l_l_i-1]);
	    --l_l_i;
	} else if ( l[l_l_i-1] == l[l_r_i+1] ) {
	    total = ls_add(total, l[l_r_i+1]);
	    total = ls_add(total, l[l_l_i-1]);
	    --l_l_i;
	    ++l_r_i;
	} else if ( l[l_l_i-1] > l[l_r_i+1] ) {
	    total = ls_add(total, l[l_l_i-1]);
	    --l_l_i;
	} else {
	    total = ls_add(total,l[l_r_i+1]);
	    ++l_r_i;
	}
    }
    CHR_POS abs_l_l_95 = start_l + l_l_i,
	abs_l_r_95 = start_l + l_r_i;

    total = r[r_max_i];
    CHR_POS r_l_i = r_max_i,
	r_r_i = r_max_i;
    t_last = bp->interval_r.i.end - bp->interval_r.i.start;

    while ( ((r_l_i > 0) || (r_r_i < t_last)) && (total < p_95) ){
	if ( r_l_i == 0 ) {
	    total = ls_add(total, r[r_r_i+1]);
	    ++r_r_i;
	} else if ( r_r_i == t_last ) {
	    total = ls_add(total, r[r_l_i-1]);
	    --r_l_i;
	} else if ( r[r_l_i-1] == r[r_r_i+1] ) {
	    total = ls_add(total, r[r_r_i+1]);
	    total = ls_add(total, r[r_l_i-1]);
	    --r_l_i;
	    ++r_r_i;
	} else if ( r[r_l_i-1] > r[r_r_i+1] ) {
	    total = ls_add(total, r[r_l_i-1]);
	    --r_l_i;
	} else {
	    total = ls_add(total,r[r_r_i+1]);
	    ++r_r_i;
	}
    }

    CHR_POS abs_r_l_95 = start_r + r_l_i,
	abs_r_r_95 = start_r + r_r_i;

    set_info("95CI", to_string(abs_l_l_95));

    // cout << "95:" <<
    // 	bp->interval_l.i.chr << ":" << abs_l_l_95 << "-"<< abs_l_r_95 <<";"<<
    // 	bp->interval_r.i.chr << ":" << abs_r_l_95 << "-"<< abs_r_r_95;

    CHR_POS open_l_start = 0,
	open_r_start = 0;

    if (bp->interval_l.i.start > 0)
	open_l_start = bp->interval_l.i.start - 1;

    if (bp->interval_r.i.start > 0)
	open_r_start = bp->interval_r.i.start - 1;

    // cout <<
    // 	interval_l.i.chr << sep <<
    // 	open_l_start << sep <<
    // 	interval_l.i.end  << sep <<
    // 	interval_r.i.chr << sep <<
    // 	open_r_start << sep<<
    // 	interval_r.i.end << sep;
    

    // cout <<
    // 	(score_l+score_r) << "\t" <<
    // 	bp->interval_l.i.strand << "\t" <<
    // 	bp->interval_r.i.strand << "\t";
    
    free(l);
    free(r);
}

void
SV_VcfVariant::
add_sample(string sample_name)
{
    if (find(samples.begin(),
	     samples.end(),
	     sample_name)
	== samples.end())
	samples.push_back(sample_name);
}

void
SV_VcfVariant::
set_info(string id)
{
    if (info != "")
	info.append(";");
    info.append(id);
}

void
SV_VcfVariant::
set_info(string id,
	 string val)
{
    if (info != "")
	info.append(";");
    info.append(id);
    info.append("=");
    info.append(val);
}

void
SV_VcfVariant::
set_sample_field(string sample,
		 string format,
		 string value)
{
    // add to format to active_formats vector
    if (find(active_formats.begin(),
	     active_formats.end(),
	     format)
	== active_formats.end())
	active_formats.push_back(format);

    // set the sample's format value
    sample_vals[sample][format] = value;
}

string
SV_VcfVariant::
get_sample_field(string sample,
		 string format)
{
    return sample_vals[sample][format];
}

string
SV_VcfVariant::
get_format_string()
{
    ostringstream fmt_ss;
    vector<string>::iterator fmt_it;
    for (fmt_it = active_formats.begin();
	 fmt_it != active_formats.end();
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
    for (fmt_it = active_formats.begin();
	 fmt_it != active_formats.end();
	 ++fmt_it) {
	// add delimiter if not first field
	if (samp_ss.str() != "")
	    samp_ss << ":";
	// store "." if no value
	if (sample_vals[sample][*fmt_it] == "")
	    samp_ss << ".";
	else
	    samp_ss << sample_vals[sample][*fmt_it];
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
    for (samp_it = samples.begin();
	 samp_it != samples.end();
	 ++samp_it)
	cout << sep <<
	    get_sample_string(*samp_it);
    
    cout << endl;
}

void
SV_VcfVariant::
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
	"##INFO=<ID=SU,Number=.,Type=Integer,Description=\"Number of pieces of evidence supporting the variant across all samples\">" << endl <<
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

    vector<string>::iterator samp_it;
    for (samp_it = samples.begin();
	 samp_it != samples.end();
	 ++samp_it) {
	cout << sep <<
	    *samp_it;
    }
    cout << endl;
}
