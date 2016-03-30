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
~SV_VcfVariant()
{
}

SV_VcfVariant::
SV_VcfVariant()
{
    is_multiline = false;
    
    chrom[LINE1] = ".";
    pos[LINE1] = 0;
    id[LINE1] = ".";
    ref[LINE1] = "N";
    alt[LINE1] = ".";
    qual[LINE1] = ".";
    filter[LINE1] = ".";
    info[LINE1] = "";
}

SV_VcfVariant::
SV_VcfVariant(SV_BreakPoint *bp,
	      int bp_id,
	      int print_prob)
{
    is_multiline = false;
    
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

    // CHROM
    chrom[LINE1] = bp->interval_l.i.chr;
    chrom[LINE2] = bp->interval_r.i.chr;
    ref[LINE1] = "N";
    ref[LINE2] = "N";
    qual[LINE1] = ".";
    qual[LINE2] = ".";
    filter[LINE1] = ".";
    filter[LINE2] = ".";
    
    // SV type
    map<string,int> uniq_strands = get_strands(bp);
    if (bp->interval_l.i.chr.compare(bp->interval_r.i.chr) != 0) {
	is_multiline = true;
	set_info(LINE1, "SVTYPE", "BND");
	set_info(LINE2, "SVTYPE", "BND");
    }
    else if (bp->type == bp->DELETION ) {
	alt[LINE1] = "<DEL>";
	set_info(LINE1, "SVTYPE", "DEL");
    }
    else if (bp->type == bp->DUPLICATION) {
	alt[LINE1] = "<DUP>";
	set_info(LINE1, "SVTYPE", "DUP");
    }
    else if (bp->type == bp->INVERSION
	     && uniq_strands.find("++") != uniq_strands.end()
	     && uniq_strands.find("--") != uniq_strands.end()) {
	alt[LINE1] = "<INV>";
	set_info(LINE1, "SVTYPE", "INV");
    }
    else {
	// set the alt for the first of two lines
	// in the multi-line variant
	is_multiline = true;
	set_info(LINE1, "SVTYPE", "BND");
	set_info(LINE2, "SVTYPE", "BND");
    }

    // INFO: STRANDS
    stringstream strands;
    stringstream strands_sec; // secondary line strands
    map<string,int>::iterator s_it;
    for (s_it = uniq_strands.begin();
	 s_it != uniq_strands.end();
	 ++s_it) {
    	if (s_it != uniq_strands.begin())
    	    strands <<  ",";
    	strands << s_it->first << ":" << s_it->second;

	// flip the orientations for the secondary breakend
	if (s_it->first == "+-")
	    strands_sec << "-+" << ":" << s_it->second;
	else if (s_it->first == "-+")
	    strands_sec << "+-" << ":" << s_it->second;
	else
	    strands_sec << s_it->first << ":" << s_it->second;
    }
    set_info(LINE1, "STRANDS", strands.str());
    if (is_multiline)
	set_info(LINE2, "STRANDS", strands_sec.str());

    if (is_multiline) {
	// for multi-line variants (non-symbolic) we need to bump
	// the position of '-' strand breakend (VCF 4.2 spec,
	// section 5.4)
	// This is because for symbolic variants, "END" refers to
	// the last deleted base in the SV, but for BND it must
	// be the immediate flanking base outside the SV.
	//
	// Symbolic format (<DEL>, <DUP>, <INV>:
	//
	//        POS      END
	//         |        |
	// =========.........=========
	//
	//
	// BND format:
	//
	//        POS       ALT
	//         |         |
	// =========.........=========
	//
	// LUMPY always reports the base to the left of the
	// breakpoint, so only need to move the '-' strand
	// positions one base to the right
	
	// move the negative strand positions to the right
	if (uniq_strands.find("--") != uniq_strands.end()) {
	    pos[LINE1] = abs_max_l + 1;
	    pos[LINE2] = abs_max_r + 1;
	}
	else if (uniq_strands.find("-+") != uniq_strands.end()) {
	    pos[LINE1] = abs_max_l + 1;
	    pos[LINE2] = abs_max_r;
	}
	else if (uniq_strands.find("++") != uniq_strands.end()) {
	    pos[LINE1] = abs_max_l;
	    pos[LINE2] = abs_max_r;
	}
	else if (uniq_strands.find("+-") != uniq_strands.end()) {
	    pos[LINE1] = abs_max_l;
	    pos[LINE2] = abs_max_r + 1;
	}
	
	set_info(LINE2, "SECONDARY");
	id[LINE1] = to_string(bp_id).append("_1");
	id[LINE2] = to_string(bp_id).append("_2");

	set_info(LINE1, "EVENT", to_string(bp_id));
	set_info(LINE2, "EVENT", to_string(bp_id));

	set_info(LINE1, "MATEID", to_string(bp_id).append("_2"));
	set_info(LINE2, "MATEID", to_string(bp_id).append("_1"));

	// ALT
	vector<string> alt_v;
	alt_v = get_alt(chrom[LINE1],
			pos[LINE1],
			ref[LINE1],
			chrom[LINE2],
			pos[LINE2],
			ref[LINE2],
			uniq_strands.begin()->first);
	
	alt[LINE1] = alt_v[LINE1];
	alt[LINE2] = alt_v[LINE2];
    }
    else {
	// POS
	pos[LINE1] = abs_max_l;

	// ID
	id[LINE1] = to_string(bp_id);
	
        // INFO: SVLEN
	int svlen;
	if (bp->type == bp->DELETION)
	    svlen = abs_max_l - abs_max_r;
	else
	    svlen = abs_max_r - abs_max_l;
	set_info(LINE1, "SVLEN", to_string(svlen));
	
	// INFO: END
	set_info(LINE1, "END", to_string(abs_max_r));
    }
        
    // INFO: CIPOS
    int cipos_l = bp->interval_l.i.start - abs_max_l;
    int cipos_r = bp->interval_l.i.end - abs_max_l;
    string cipos = to_string(cipos_l);
    cipos.append(",");
    cipos.append(to_string(cipos_r));

    int ciend_l = bp->interval_r.i.start - abs_max_r;
    int ciend_r = bp->interval_r.i.end - abs_max_r;
    string ciend = to_string(ciend_l);
    ciend.append(",");
    ciend.append(to_string(ciend_r));

    set_info(LINE1, "CIPOS", cipos);
    set_info(LINE1, "CIEND", ciend);
   
    if (is_multiline) {
	set_info(LINE2, "CIPOS", ciend);
	set_info(LINE2, "CIEND", cipos);
    }

    // INFO: Get the area that includes 95% of the probabitliy
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

    // CHR_POS abs_l_l_95 = start_l + l_l_i,
    // abs_l_r_95 = start_l + l_r_i;
    int cipos95_l = start_l + l_l_i - abs_max_l;
    int cipos95_r = start_l + l_r_i - abs_max_l;
    string cipos95 = to_string(cipos95_l);
    cipos95.append(",");
    cipos95.append(to_string(cipos95_r));

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

    // CHR_POS abs_r_l_95 = start_r + r_l_i,
    // 	abs_r_r_95 = start_r + r_r_i;
    int ciend95_l = start_r + r_l_i - abs_max_r;
    int ciend95_r = start_r + r_r_i - abs_max_r;
    string ciend95 = to_string(ciend95_l);
    ciend95.append(",");
    ciend95.append(to_string(ciend95_r));

    set_info(LINE1, "CIPOS95", cipos95);
    set_info(LINE1, "CIEND95", ciend95);
    
    if (is_multiline) {
	set_info(LINE2, "CIPOS95", ciend95);
	set_info(LINE2, "CIEND95", cipos95);
    }

    // INFO: IMPRECISE (based on 95% confidence interval)
    if (cipos95_l != 0 ||
    	cipos95_r != 0 ||
    	ciend95_l != 0 ||
    	ciend95_r != 0) {
    	set_info(LINE1, "IMPRECISE");
    	if (is_multiline)
    	    set_info(LINE2, "IMPRECISE");
    }

    // FORMAT: initialize appropriate sample fields
    vector<string> ev_v;
    ev_v.push_back("PE");
    ev_v.push_back("SR");
    ev_v.push_back("BD");

    vector<string>::iterator samp_it;
    for (samp_it = samples.begin();
	 samp_it != samples.end();
	 ++samp_it) {
	set_sample_field(*samp_it, "GT", "./.");
	set_sample_field(*samp_it, "SU", "0");
	for (vector<string>::iterator it = ev_v.begin();
	     it != ev_v.end();
	     ++it) {
	    if (find(active_formats.begin(),
		     active_formats.end(),
		     *it)
		!= active_formats.end())
		set_sample_field(*samp_it, *it, "0");
	}
    }

    // FORMAT: update variant support with evidence counts
    map<string,int> var_supp;
    map<int, int>::iterator ids_it;
    for (ids_it = bp->ev_ids.begin();
	 ids_it != bp->ev_ids.end();
	 ++ids_it) {
	string samp = SV_EvidenceReader::sample_names[ids_it->first];
	string ev_type = SV_EvidenceReader::ev_types[ids_it->first];

	int samp_supp = atoi(get_sample_field(samp, "SU").c_str());
	int ev = atoi(get_sample_field(samp, ev_type).c_str());
	int new_ev = ids_it->second;

	ev += new_ev;
	samp_supp += new_ev;

	var_supp["SU"] += new_ev;
	var_supp[ev_type] += new_ev;

	set_sample_field(samp, ev_type, to_string(ev));
	set_sample_field(samp, "SU", to_string(samp_supp));
    }
    
    // INFO: add variant support across all samples
    set_info(LINE1, "SU", to_string(var_supp["SU"]));
    if (is_multiline)
	set_info(LINE2, "SU", to_string(var_supp["SU"]));
    for (vector<string>::iterator ev_it = ev_v.begin();
	 ev_it != ev_v.end();
	 ++ev_it) {
	if (find(active_formats.begin(),
		 active_formats.end(),
		 *ev_it)
	    != active_formats.end()) {
	    set_info(LINE1, *ev_it, to_string(var_supp[*ev_it]));
	    if (is_multiline)
		set_info(LINE2, *ev_it, to_string(var_supp[*ev_it]));
	}
    }

    // INFO: PRPOS, PREND
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
	set_info(LINE1, "PRPOS", lp);
	set_info(LINE1, "PREND", rp);

	if (is_multiline) {
	    set_info(LINE2, "PRPOS", rp);
	    set_info(LINE2, "PREND", lp);
	}
    }

    free(l);
    free(r);
}

map<string,int>
SV_VcfVariant::
get_strands(SV_BreakPoint *bp)
{
    vector<SV_Evidence*>::iterator it;
    vector<SV_BreakPoint *> bps;
    map<string,int> uniq_strands;
    for (it = bp->evidence.begin();
	 it < bp->evidence.end();
	 ++it) {
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

    vector<SV_BreakPoint *>::iterator bp_it;
    for (bp_it = bps.begin(); bp_it < bps.end(); ++bp_it) {
	SV_BreakPoint *tmp_bp = *bp_it;
	delete tmp_bp;
    }
    return uniq_strands;
}

vector<string>
SV_VcfVariant::
get_alt(string l_chrom,
	CHR_POS l_pos,
	string l_ref,
	string r_chrom,
	CHR_POS r_pos,
	string r_ref,
	string strands)
{
    vector<string> alt;
    stringstream l_alt, r_alt;

    if (strands == "+-") {
	l_alt << l_ref <<
	    "[" << r_chrom <<
	    ":" << r_pos << "[";
	r_alt << "]" <<
	    l_chrom << ":" <<
	    l_pos << "]" << r_ref;
    }
    else if (strands == "-+") {
	l_alt << "]" <<
	    r_chrom << ":" <<
	    r_pos << "]" << l_ref;
	r_alt << r_ref <<
	    "[" << l_chrom <<
	    ":" << l_pos << "[";
    }
    else if (strands == "++") {
	l_alt << l_ref <<
	    "]" << r_chrom <<
	    ":" << r_pos << "]";
	r_alt << r_ref <<
	    "]" << l_chrom <<
	    ":" << l_pos << "]";
    }
    else if (strands == "--") {
	l_alt << "[" <<
	    r_chrom << ":" <<
	    r_pos << "[" << l_ref;
	r_alt << "[" <<
	    l_chrom << ":" <<
	    l_pos << "[" << r_ref;
    }

    alt.push_back(l_alt.str());
    alt.push_back(r_alt.str());
    
    return alt;
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
set_info(int line,
	 string id)
{
    if (info[line] != "")
	info[line].append(";");
    info[line].append(id);
}

void
SV_VcfVariant::
set_info(int line,
	 string id,
	 string val)
{
    if (info[line] != "")
	info[line].append(";");
    info[line].append(id);
    info[line].append("=");
    info[line].append(val);
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
    int num_lines = is_multiline + 1;

    for (int i = 0; i < num_lines; ++i) {
	cout << chrom[i] << sep <<
	    pos[i] << sep <<
	    id[i] << sep <<
	    ref[i] << sep <<
	    alt[i] << sep <<
	    qual[i] << sep <<
	    filter[i] << sep <<
	    info[i] << sep <<
	    format;
	vector<string>::iterator samp_it;
	for (samp_it = samples.begin();
	     samp_it != samples.end();
	     ++samp_it)
	    cout << sep <<
		get_sample_string(*samp_it);
    
	cout << endl;
    }
}

void
SV_VcfVariant::
print_header()
{
    string sep = "\t";

    cout <<
	"##fileformat=VCFv4.2" << endl <<
	"##source=LUMPY" << endl <<
	"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl <<
	"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl <<
	"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl <<
	"##INFO=<ID=STRANDS,Number=.,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">" << endl <<
	"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << endl <<
	"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << endl <<
	"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << endl <<
	"##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description=\"Confidence interval (95%) around POS for imprecise variants\">" << endl <<
	"##INFO=<ID=CIEND95,Number=2,Type=Integer,Description=\"Confidence interval (95%) around END for imprecise variants\">" << endl <<
	"##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">" << endl <<
	"##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">" << endl <<
	"##INFO=<ID=SECONDARY,Number=0,Type=Flag,Description=\"Secondary breakend in a multi-line variants\">" << endl <<
	"##INFO=<ID=SU,Number=.,Type=Integer,Description=\"Number of pieces of evidence supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=PE,Number=.,Type=Integer,Description=\"Number of paired-end reads supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=SR,Number=.,Type=Integer,Description=\"Number of split reads supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=BD,Number=.,Type=Integer,Description=\"Amount of BED evidence supporting the variant across all samples\">" << endl <<
	"##INFO=<ID=EV,Number=.,Type=String,Description=\"Type of LUMPY evidence contributing to the variant call\">" << endl <<
	"##INFO=<ID=PRPOS,Number=.,Type=String,Description=\"LUMPY probability curve of the POS breakend\">" << endl <<
	"##INFO=<ID=PREND,Number=.,Type=String,Description=\"LUMPY probability curve of the END breakend\">" << endl <<
	"##ALT=<ID=DEL,Description=\"Deletion\">" << endl <<
	"##ALT=<ID=DUP,Description=\"Duplication\">" << endl <<
	"##ALT=<ID=INV,Description=\"Inversion\">" << endl <<
	"##ALT=<ID=DUP:TANDEM,Description=\"Tandem duplication\">" << endl <<
	"##ALT=<ID=INS,Description=\"Insertion of novel sequence\">" << endl <<
	"##ALT=<ID=CNV,Description=\"Copy number variable region\">" << endl <<
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl <<
	"##FORMAT=<ID=SU,Number=1,Type=Integer,Description=\"Number of pieces of evidence supporting the variant\">" << endl <<
	"##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Number of paired-end reads supporting the variant\">" << endl <<
	"##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of split reads supporting the variant\">" << endl <<
    	"##FORMAT=<ID=BD,Number=1,Type=Integer,Description=\"Amount of BED evidence supporting the variant\">" << endl;

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

template <typename Type>
string
SV_VcfVariant::
to_string(Type t)
{
  stringstream convert;
  convert << t;
  return convert.str();
}
