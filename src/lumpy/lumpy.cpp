/*****************************************************************************
 * lumpy.cpp
 * (c) 2012 - Ryan M. Layer
 * Hall Laboratory
 * Quinlan Laboratory
 * Department of Computer Science
 * Department of Biochemistry and Molecular Genetics
 * Department of Public Health Sciences and Center for Public Health Genomics,
 * University of Virginia
 * rl6sf@virginia.edu
 *
 * Licenced under the GNU General Public License 2.0 license.
 * ****************************************************************************/

//{{{ includes
#include "version.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BamAncillary.h"
#include "bedFile.h"
#include "bedFilePE.h"
#include "sequenceUtils.h"
#include "SV_Evidence.h"
#include "SV_Pair.h"
#include "SV_BreakPoint.h"
#include "SV_Bedpe.h"
#include "SV_SplitRead.h"
#include "SV_SplitReadReader.h"
#include "SV_PairReader.h"
#include "SV_BedpeReader.h"
#include "SV_Tools.h"

#include "ucsc_bins.hpp"
#include "log_space.h"

using namespace BamTools;

#include <exception>
#include <vector>
#include <map>
#include <algorithm>  
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <limits.h>
#include <math.h>

#include <gsl_statistics_int.h>


using namespace std;


// define our program name
#define PROGRAM_NAME "**********"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)
//}}}

//{{{ forward declarations
void ShowHelp(void);

static inline int strnum_cmp(const char *a, const char *b);
//}}}

//{{{ void ShowHelp(void)
void ShowHelp(void)
{
    cerr << endl << "Program: " << PROGRAM_NAME << " (v 0.1)" <<
			endl <<
		"Author:  Ryan Layer (rl6sf@virginia.edu)" << endl <<
		"Summary: Find structural variations in various signals." << endl <<
			endl <<
		"Usage:   " << PROGRAM_NAME << " [OPTIONS] " << endl << endl <<
		"Options: " << endl <<
		"\t-e"		"\tShow evidnece for each call" << endl <<
		"\t-mw"		"\tminimum weight for a call" << endl <<
		"\t-tt"		"\ttrim threshold" << endl <<
		"\t-sr"		"\tbam_file:<file name>," << endl <<
					"\t\tback_distance:<distance>" << endl << 
					"\t\tmin_mapping_threshold:<mapping quality>" << endl << 
					"\t\tweight:<sample weight>" << endl << 
					"\t\tid:<sample id>" << endl << 
					endl <<
		"\t-pe"		"\tbam_file:<file name>," << endl <<
					"\t\thisto_file:<file name>," << endl <<
					"\t\tmean:<value>," << endl << 
					"\t\tstdev:<value>," << endl <<
					"\t\tread_length:<length>," << endl <<
					"\t\tmin_non_overlap:<length>," << endl <<
					"\t\tdiscordant_z:<z value>," << endl << 
					"\t\tback_distance:<distance>" << endl <<
					"\t\tmin_mapping_threshold:<mapping quality>" << endl << 
					"\t\tweight:<sample weight>" << endl << 
					"\t\tid:<sample id>" << endl << 
					endl << 
		"\t-pe"		"\tbedpe_file:<bedpe file>," << endl <<
					"\t\tdistro_file:<distro_file>," << endl <<
					"\t\tback_distance:<back distance>" << endl <<
					"\t\tweight:<sample weight>" << endl << 
					"\t\tid:<sample id>" << endl << 
					endl;
    // end the program here
    exit(1);
}
//}}}

//{{{ static inline int strnum_cmp(const char *a, const char *b)
//read name str cmp
static inline int strnum_cmp(const char *a, const char *b)
{
	char *pa, *pb;
	pa = (char*)a; pb = (char*)b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			long ai, bi;
			ai = strtol(pa, &pa, 10);
			bi = strtol(pb, &pb, 10);
			if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
		} else {
			if (*pa != *pb) break;
			++pa; ++pb;
		}
	}
	if (*pa == *pb)
		return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
	return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}
//}}}

//{{{ void parse_cmd_line_args(int *i,
void
parse_cmd_line_args(int *i,
					int argc,
					char **argv,
					SV_EvidenceReader *e_r,
					vector<SV_EvidenceReader*> &evidence_readers,
					string e_name)
{
	if ((*i+1) < argc) {
		char *params = argv[*i + 1];
		char *param_val, *brka, *brkb;

		for (	param_val = strtok_r(params, ",", &brka);
				param_val;
				param_val = strtok_r(NULL, ",", &brka)) {   
			char *param = strtok_r(param_val, ":", &brkb);
			char *val = strtok_r(NULL, ":", &brkb);

			if (val == NULL) {
				cerr << "Parameter requied for " << param << endl;
				ShowHelp();
			}

			if ( ! e_r->add_param(param, val) ) {
				cerr << "Unknown pair end parameter:" << param << endl;
				ShowHelp();
			}
		}
	}

	string msg = e_r->check_params();
	if ( msg.compare("") == 0 ) {
		evidence_readers.push_back(e_r);
	} else {
		cerr << "missing " << e_name << " parameters:" << msg << endl;
		ShowHelp();
	}

	i++;

}
//}}}

//{{{ int main(int argc, char* argv[]) {
int main(int argc, char* argv[])
{

	//{{{ setup
	double trim_threshold = 1e-10;
	double merge_threshold = 1e-10;
	int min_weight = 0;
	bool show_evidence = false;
	//}}}

    //{{{ check to see if we should print out some help
	if (argc == 1) 
		ShowHelp();

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
           (PARAMETER_CHECK("--help", 6, parameterLength))) {
			ShowHelp();
        }
    }
	//}}}

	//{{{ do some parsing (all of these parameters require 2 strings)
	vector<SV_EvidenceReader*> evidence_readers;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

		if(PARAMETER_CHECK("-pe", 3, parameterLength)) {
			//{{{

			SV_PairReader *pe_r = new SV_PairReader();

			/*
			parse_cmd_line_args(&i,
								argc,
								argv,
								pe_r,
								evidence_readers,
								"pair-end");
			*/
			if ((i+1) < argc) {
				char *params = argv[i + 1];
				char *param_val, *brka, *brkb;

				for (	param_val = strtok_r(params, ",", &brka);
						param_val;
						param_val = strtok_r(NULL, ",", &brka)) {   
					char *param = strtok_r(param_val, ":", &brkb);
					char *val = strtok_r(NULL, ":", &brkb);

					if (val == NULL) {
						cerr << "Parameter requied for " << param << endl;
						ShowHelp();
					}

					if ( ! pe_r->add_param(param, val) ) {
						cerr << "Unknown pair end parameter:" << param << endl;
						ShowHelp();
					}
				}
			}

			string msg = pe_r->check_params();
			if ( msg.compare("") == 0 ) {
				// Add to list of readers
				evidence_readers.push_back(pe_r);
				// Set the ditro map
				pe_r->set_statics();
				SV_Evidence::distros[pe_r->sample_id] = 
					pair<log_space*,log_space*>(
							SV_Pair::get_bp_interval_probability('+'),
							SV_Pair::get_bp_interval_probability('-'));
			} else {
				cerr << "missing pair end parameters:" << msg << endl;
				ShowHelp();
			}

			i++;
			//}}}
		}
	

		else if(PARAMETER_CHECK("-bedpe", 6, parameterLength)) {
			//{{{
			
			SV_BedpeReader *be_r = new SV_BedpeReader();

			if ((i+1) < argc) {
				char *params = argv[i + 1];
				char *param_val, *brka, *brkb;

				for (	param_val = strtok_r(params, ",", &brka);
						param_val;
						param_val = strtok_r(NULL, ",", &brka)) {   
					char *param = strtok_r(param_val, ":", &brkb);
					char *val = strtok_r(NULL, ":", &brkb);

					if (val == NULL) {
						cerr << "Parameter requied for " << param << endl;
						ShowHelp();
					}
	
					if ( ! be_r->add_param(param, val) ) {
						cerr << "Unknown bedpe parameter:" << param << endl;
						ShowHelp();
					}
				}
			}

			string msg = be_r->check_params();
			if ( msg.compare("") == 0 ) {
				evidence_readers.push_back(be_r);
				SV_Evidence::distros[be_r->sample_id] = 
					pair<log_space*,log_space*>(
							SV_Bedpe::get_bp_interval_probability('+'),
							SV_Bedpe::get_bp_interval_probability('-'));

			} else {
				cerr << "missing pair end parameters:" << msg << endl;
				ShowHelp();
			}

			i++;
			//}}}
		}

		else if(PARAMETER_CHECK("-sr", 3, parameterLength)) {
			//{{{
			SV_SplitReadReader *sr_r = new SV_SplitReadReader();

			if ((i+1) < argc) {
				char *params = argv[i + 1];
				char *param_val, *brka, *brkb;

				for (	param_val = strtok_r(params, ",", &brka);
						param_val;
						param_val = strtok_r(NULL, ",", &brka)) {   
					char *param = strtok_r(param_val, ":", &brkb);
					char *val = strtok_r(NULL, ":", &brkb);

					if (val == NULL) {
						cerr << "Parameter requied for " << param << endl;
						ShowHelp();
					}
	
					if ( ! sr_r->add_param(param, val) ) {
						cerr << "Unknown pair end parameter:" << param << endl;
						ShowHelp();
					}
				}
			}

			string msg = sr_r->check_params();
			if ( msg.compare("") == 0 ) {
				evidence_readers.push_back(sr_r);
				SV_Evidence::distros[sr_r->sample_id] = 
					pair<log_space*,log_space*>(
							SV_SplitRead::get_bp_interval_probability('+'),
							SV_SplitRead::get_bp_interval_probability('-'));

			} else {
				cerr << "missing pair end parameters:" << msg << endl;
				ShowHelp();
			}

			i++;
			//}}}
		}

        else if(PARAMETER_CHECK("-tt", 3, parameterLength)) {
            if ((i+1) < argc) {
                trim_threshold = atof(argv[i + 1]);
                i++;
			}
		}

        else if(PARAMETER_CHECK("-mt", 3, parameterLength)) {
            if ((i+1) < argc) {
                merge_threshold = atof(argv[i + 1]);
                i++;
			}
		}
        else if(PARAMETER_CHECK("-mw", 3, parameterLength)) {
            if ((i+1) < argc) {
                min_weight = atoi(argv[i + 1]);
                i++;
			}
		}

        else if(PARAMETER_CHECK("-e", 2, parameterLength)) {
			show_evidence = true;
		}

		/*
        else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
            if ((i+1) < argc) {
				store_evidence = true;
                evidence_db = argv[i + 1];
				i++;
			}
		}
		*/

        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " <<
					argv[i] << " *****" << endl << endl;
			ShowHelp();
        }
	}

	if (min_weight == 0) {
		cerr << endl << "*****ERROR: must set min weight *****" << endl << endl;
		ShowHelp();
	}


	//}}}

	UCSCBins<SV_BreakPoint*> r_bin;

	SV_BreakPoint::p_trim_threshold = trim_threshold;
	SV_BreakPoint::p_merge_threshold = merge_threshold;

	vector<SV_EvidenceReader*>::iterator i_er;

	//{{{ initialize all input files
	for (	i_er = evidence_readers.begin();
			i_er != evidence_readers.end();
			++i_er) {
		SV_EvidenceReader *er = *i_er;
		er->initialize();
	}
	//}}}
	
	bool has_next = false;

	//{{{ Test if there lines to process in each input file
	for (	i_er = evidence_readers.begin();
			i_er != evidence_readers.end();
			++i_er) {
		SV_EvidenceReader *er = *i_er;
		has_next = has_next || er->has_next();
	}
	//}}}

	while ( has_next ) {
		string min_chr = "";
		//{{{ find min_chr among all input files
		for (	i_er = evidence_readers.begin();
				i_er != evidence_readers.end();
				++i_er) {
			SV_EvidenceReader *er = *i_er;

			if ( er->has_next() ) {
				string curr_chr = er->get_curr_chr();

				if ( ( min_chr.compare("") == 0 ) ||
					 ( curr_chr.compare(min_chr) < 0 ) ) {
					min_chr = curr_chr;
				}
			}
		}
		//}}}

		//cerr << min_chr << endl;

		// {{{ loop through each input file and process min_chr
		for (	i_er = evidence_readers.begin();
				i_er != evidence_readers.end();
				++i_er) {

			SV_EvidenceReader *er = *i_er;

			if ( er->has_next() ) {
				string curr_chr = er->get_curr_chr();
				if ( curr_chr.compare(min_chr) <= 0 )  {
					er->set_statics();
					er->process_input_chr(curr_chr, r_bin);
				}
			}
		}
		//}}}

		//{{{ get breakpoints where both sides are in this chrm

		vector< UCSCElement<SV_BreakPoint*> > values = r_bin.values(min_chr);
		vector< UCSCElement<SV_BreakPoint*> >::iterator it;

#if 0
		//{{{ sort and uniq
		cerr << "Sorting..." << endl;
		sort(values.begin(), values.end(),
				UCSCElement<SV_BreakPoint*>::sort_ucscelement_by_value);
		cerr << "Done." << endl;
		cerr << "Unique..." << endl;
		it = unique(values.begin(),
					values.end(),
					UCSCElement<SV_BreakPoint*>::compare_ucscelement_by_value);
		cerr << "Done." << endl;

		values.resize( it - values.begin() );
		//}}}
#endif
		
		vector< UCSCElement<SV_BreakPoint*> > to_remove;
		//{{{ get breakpoints
		for (it = values.begin(); it < values.end(); ++it) {
			SV_BreakPoint *bp = it->value;

			vector< UCSCElement<SV_BreakPoint*> > gets = 
					r_bin.get(it->chr,
							  it->start,
							  it->end,
							  it->strand,
							  false);

			if ( ( bp->interval_l.i.chr.compare(min_chr) == 0 ) &&
				 ( bp->interval_r.i.chr.compare(min_chr) == 0 ) ) {

				if ( bp->weight >= min_weight ) {
					 
					bp->trim_intervals();
					bp->print_bedpe(-1);
					if (show_evidence)
						bp->print_evidence("\t");
				}

				r_bin.remove(*it, false, false, true);
				bp->free_evidence();
				delete bp;
			}
		}
		//}}}

		//}}}

		has_next = false;
		//{{{ Test if there is still input lines
		for (	i_er = evidence_readers.begin();
				i_er != evidence_readers.end();
				++i_er) {
			SV_EvidenceReader *er = *i_er;
			has_next = has_next || er->has_next();
		}
		//}}}
	}
	
	//{{{ terminate input files
	for (	i_er = evidence_readers.begin();
			i_er != evidence_readers.end();
			++i_er) {
		SV_EvidenceReader *er = *i_er;
		er->terminate();
	}
	//}}}
	
	//cerr << "Done with input files" << endl;

	//{{{ Call break points
	vector< UCSCElement<SV_BreakPoint*> > values = r_bin.values();
	vector< UCSCElement<SV_BreakPoint*> >::iterator it;

		for (it = values.begin(); it != values.end(); ++it) {
			SV_BreakPoint *bp = it->value;

			if ( bp->weight >= min_weight ) {
				 
				bp->trim_intervals();
				bp->print_bedpe(-1);
				if (show_evidence)
					bp->print_evidence("\t");
			}

			bp->free_evidence();
			delete bp;
		}

	//}}}
	return 1;
}
