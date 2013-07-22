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
#include "api/BamMultiReader.h"
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
#include "SV_BamReader.h"
#include "SV_InterChromBamReader.h"
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
#include <cstdio>
#include <cstdlib>
#include <time.h>


#include <gsl_statistics_int.h>


using namespace std;


// define our program name
#define PROGRAM_NAME "**********"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)
//}}}

//{{{ forward declarations
void ShowHelp(void);
//}}}

//{{{ void ShowHelp(void)
void ShowHelp(void)
{
    cerr << endl << "Program: " << PROGRAM_NAME << " (v 0.1.1)" <<
                endl <<
    "Author:  Ryan Layer (rl6sf@virginia.edu)" << endl <<
    "Summary: Find structural variations in various signals." << endl <<
            endl <<
    "Usage:   " << PROGRAM_NAME << " [OPTIONS] " << endl << endl <<
    "Options: " << endl <<
        "\t-e"	"\tShow evidnece for each call" << endl <<
        "\t-w"	"\tFile read windows size (default 1000000)" << endl <<
        "\t-mw"	"\tminimum weight for a call" << endl <<
        "\t-tt"	"\ttrim threshold" << endl <<
        "\t-x"	"\texclude file bed file" <<  endl <<
        "\t-t"	"\ttemp file prefix, must be to a writeable directory" << 
    endl <<
        "\t-sr"     "\tbam_file:<file name>," << endl <<
            "\t\tback_distance:<distance>" << endl << 
            "\t\tmin_mapping_threshold:<mapping quality>" << endl << 
            "\t\tweight:<sample weight>" << endl << 
            "\t\tid:<sample id>" << endl << 
            "\t\tmin_clip:<minimum clip length>" << endl << 
    endl <<
        "\t-pe"     "\tbam_file:<file name>," << endl <<
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
        "\t-bedpe"  "\tbedpe_file:<bedpe file>," << endl <<
            "\t\tdistro_file:<distro_file>," << endl <<
            "\t\tback_distance:<back distance>" << endl <<
            "\t\tweight:<sample weight>" << endl << 
            "\t\tid:<sample id>" << endl << 
    endl;
    // end the program here
    exit(1);
}
//}}}

int main(int argc, char* argv[])
{

    //{{{ setup
    double trim_threshold = 1e-10;
    double merge_threshold = 1e-10;
    int min_weight = 0;
    bool show_evidence = false;
    bool has_bams = false, has_bedpes = false;
    CHR_POS window_size = 1000000;
    string inter_chrom_file_prefix = "./";
    int call_id = 0;
    bool has_next = false;
    vector<SV_EvidenceReader*>::iterator i_er;
    UCSCBins<SV_BreakPoint*> r_bin;
    srand (time(NULL));
    string exclude_bed_file;
    bool has_exclude = false;
    //vector<string> bam_files;
    //}}}

    //{{{ check to see if we should print out some help
    if (argc == 1) 
        ShowHelp();

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if( (PARAMETER_CHECK("-h", 2, parameterLength)) ||
            (PARAMETER_CHECK("--help", 6, parameterLength))) {
            ShowHelp();
        }
    }
	//}}}

    //{{{ do some parsing and setup
    vector<SV_EvidenceReader*> evidence_readers;
    map<string, SV_EvidenceReader*> bam_evidence_readers;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-pe", 3, parameterLength)) {
                //{{{
                has_bams = true;
                SV_PairReader *pe_r = new SV_PairReader();

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
                        // Set the ditro map
                        
                        pe_r->initialize();
                        SV_Evidence::distros[pe_r->sample_id] = 
                                pair<log_space*,log_space*>(
                                        SV_Pair::get_bp_interval_probability('+',
                                                                                                                 pe_r->distro_size,
                                                                                                                 pe_r->distro),
                                        SV_Pair::get_bp_interval_probability('-',
                                                                                                                 pe_r->distro_size,
                                                                                                                 pe_r->distro));
                        SV_Evidence::distros_size[pe_r->sample_id] = pe_r->distro_size;
                } else {
                        cerr << "missing pair end parameters:" << msg << endl;
                        ShowHelp();
                }

                bam_evidence_readers[pe_r->get_source_file_name()] = pe_r;

                i++;
                //}}}
        }

        else if(PARAMETER_CHECK("-bedpe", 6, parameterLength)) {
                //{{{
                has_bedpes = true;	
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
                        be_r->initialize();
                        SV_Evidence::distros[be_r->sample_id] = 
                                pair<log_space*,log_space*>(
                                        SV_Bedpe::
                                        get_bp_interval_probability('+',
                                                                                            be_r->distro_size,
                                                                                            be_r->distro),
                                        SV_Bedpe::
                                        get_bp_interval_probability('-',
                                                                                            be_r->distro_size,
                                                                                            be_r->distro));
                        SV_Evidence::distros_size[be_r->sample_id] = be_r->distro_size;
                        evidence_readers.push_back(be_r);
                } else {
                        cerr << "missing bedpe parameters:" << msg << endl;
                        ShowHelp();
                }

                i++;
                //}}}
        }

        else if(PARAMETER_CHECK("-sr", 3, parameterLength)) {
                //{{{
                has_bams = true;
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
                        sr_r->initialize();
                        SV_Evidence::distros[sr_r->sample_id] = 
                                pair<log_space*,log_space*>(
                                                SV_SplitRead::
                                                get_bp_interval_probability('+',
                                                                                                        sr_r->back_distance),
                                                SV_SplitRead::
                                                get_bp_interval_probability('-',
                                                                                                        sr_r->back_distance));

                        SV_Evidence::distros_size[sr_r->sample_id] = 
                                        sr_r->back_distance * 2 + 1;
                } else {
                        cerr << "missing split read parameters:" << msg << endl;
                        ShowHelp();
                }

                //bam_files.push_back(sr_r->get_source_file_name());
                bam_evidence_readers[sr_r->get_source_file_name()] = sr_r;

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

        else if(PARAMETER_CHECK("-w", 2, parameterLength)) {
            if ((i+1) < argc) {
                window_size = atoi(argv[i + 1]);
                i++;
            }
        }

        else if(PARAMETER_CHECK("-x", 2, parameterLength)) {
            if ((i+1) < argc) {
                exclude_bed_file = argv[i + 1];
                has_exclude = true;
                i++;
            }
        }

        else if(PARAMETER_CHECK("-t", 2, parameterLength)) {
            if ((i+1) < argc) {
		inter_chrom_file_prefix = argv[i + 1];
                i++;
            }
        }

        else if(PARAMETER_CHECK("-e", 2, parameterLength)) {
	    show_evidence = true;
        }

        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " <<
                    argv[i] << " *****" << endl << endl;
	    ShowHelp();
        }
    }

    if (min_weight == 0) {
        cerr << endl << "*****ERROR: must set min weight *****" << 
                endl << endl;
        ShowHelp();
    }

    SV_BreakPoint::p_trim_threshold = trim_threshold;
    SV_BreakPoint::p_merge_threshold = merge_threshold;

    // append rand number to the temp inter-chrom file
    inter_chrom_file_prefix = inter_chrom_file_prefix + ToString(rand());

    if (has_exclude) 
        parse_exclude_file(exclude_bed_file, SV_Evidence::exclude_regions);

    //}}} end parsing 

    //{{{ put bams into evidence reader
    if (has_bams) {
        SV_BamReader *bam_r = new SV_BamReader(&bam_evidence_readers);
        bam_r->set_inter_chrom_file_name(inter_chrom_file_prefix + ".bam");
        bam_r->initialize();
        evidence_readers.push_back(bam_r);
    }
    //}}}

    //{{{ Test if there lines to process in each input file
    for ( i_er = evidence_readers.begin();
          i_er != evidence_readers.end();
          ++i_er) {
        SV_EvidenceReader *er = *i_er;
        has_next = has_next || er->has_next();
    }
    //}}}

    //{{{ process the intra-chrom events that were saved to a file
    CHR_POS max_pos = 0;
    string last_min_chr = "";
    while ( has_next ) {
        string min_chr = "";
        //{{{ find min_chr among all input files
        for ( i_er = evidence_readers.begin();
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

        //{{{ if the chrome switches, reset the max_pos
        if (last_min_chr.compare(min_chr) != 0) {
            max_pos = window_size;
            last_min_chr = min_chr;
        }
        //}}}

        bool input_processed = true;

        while (input_processed) {
            input_processed = false;
                            
            //{{{ read the files
            for ( i_er = evidence_readers.begin();
                  i_er != evidence_readers.end();
                  ++i_er) {

                SV_EvidenceReader *er = *i_er;

                if ( er->has_next() ) {
                    string curr_chr = er->get_curr_chr();
                    CHR_POS curr_pos = er->get_curr_pos();
                    if ( ( curr_chr.compare(min_chr) <= 0 ) &&
                         ( curr_pos < max_pos) ) {
                        er->process_input_chr_pos(curr_chr, max_pos, r_bin);
                        input_processed = true;
                    } 
                }
            }
            //}}}

            //{{{ call breakpoints
            vector< UCSCElement<SV_BreakPoint*> > values = 
                    r_bin.values(min_chr, max_pos);

            vector< UCSCElement<SV_BreakPoint*> >::iterator it;

            for (it = values.begin(); it < values.end(); ++it) {
                SV_BreakPoint *bp = it->value;

                // Make sure both ends of the bp are less than or equal to the
                // current chrom
                if ( bp->weight >= min_weight ) {
                    //bp->do_it();
                    bp->trim_intervals();
                    bp->print_bedpe(++call_id);
                    if (show_evidence)
                        bp->print_evidence("\t");
                }

                if (r_bin.remove(*it, false, false, true) != 0) {
                    cerr << "Error removing element:" << *bp << endl;
                    abort();
                }
                bp->free_evidence();
                delete bp;
            }
            //}}}

            // move the window
            max_pos = max_pos *2;
        }

        //{{{ Test if there is still input lines
        has_next = false;
        for ( i_er = evidence_readers.begin();
              i_er != evidence_readers.end();
              ++i_er) {
            SV_EvidenceReader *er = *i_er;
            has_next = has_next || er->has_next();
        }
        //}}}
    }
    //}}}
	
    //{{{ terminate input files
    for ( i_er = evidence_readers.begin();
          i_er != evidence_readers.end();
          ++i_er) {
        SV_EvidenceReader *er = *i_er;
        er->terminate();
    }
    //}}}
	
    //{{{ Call remaining intra breakpoints
    vector< UCSCElement<SV_BreakPoint*> > values = r_bin.values();
    vector< UCSCElement<SV_BreakPoint*> >::iterator it;

    for ( it = values.begin();
          it != values.end(); ++it) {
          SV_BreakPoint *bp = it->value;

        if ( bp->weight >= min_weight ) {
	    //bp->do_it();
            bp->trim_intervals();
            bp->print_bedpe(++call_id);
            if (show_evidence)
                bp->print_evidence("\t");
        }

        if (r_bin.remove(*it, false, false, true) != 0) {
	    cerr << "Error removing element" << endl;
            abort();
        }
        bp->free_evidence();
        delete bp;
    }
    //}}}
	
    //{{{ process the inter-chrom events that were saved to a file
    string intra_bam_file_name = inter_chrom_file_prefix + ".bam";
    ifstream intra_bam_file( intra_bam_file_name.c_str() );
    if (intra_bam_file.good()) {
        intra_bam_file.close();

        sort_inter_chrom_bam( inter_chrom_file_prefix + ".bam",
			      inter_chrom_file_prefix + ".sort.bam");

        SV_InterChromBamReader *ic_r = new SV_InterChromBamReader(
				inter_chrom_file_prefix + ".sort.bam",
				&bam_evidence_readers);
        ic_r->initialize();

        vector<SV_EvidenceReader*> inter_chrom_evidence_readers;
        inter_chrom_evidence_readers.push_back(ic_r);

        // There are two files containg all of the inter-chrom events, one bam
        // and one bedpe, each line in the file corresponds to the properies
        // set in one of the readers.  Each line has a "LS" (lumpy source)
        // property that gives its source file name.  Using that entry, the
        // line will be sent to the reader for processing.
		
        // get new evidence readers for both bedpe and bam inter-chrom
        // readers
        has_next = true;	

        int32_t last_min_primary_refid = -1;
        int32_t last_min_secondary_refid = -1;
        max_pos = 0;
        while ( has_next ) {
            string min_primary_chr = "";
            string min_secondary_chr = "";
            int32_t min_primary_refid = -1;
            int32_t min_secondary_refid = -1;

            //{{{ find min_chr pair among all input files
            for ( i_er = inter_chrom_evidence_readers.begin();
		  i_er != inter_chrom_evidence_readers.end();
		  ++i_er) {
		SV_EvidenceReader *er = *i_er;

		if ( er->has_next() ) {
                    int32_t curr_primary_refid = er->get_curr_primary_refid();
                    int32_t curr_secondary_refid = 
                            er->get_curr_secondary_refid();

                    if ( (( min_primary_refid == -1 ) &&
                          ( min_secondary_refid == -1 )) ||
                         (( curr_primary_refid < min_primary_refid)  && 
			  ( curr_secondary_refid < min_secondary_refid)) ){
                        min_primary_refid = curr_primary_refid;
		        min_secondary_refid = curr_secondary_refid;
		        min_secondary_chr = er->get_curr_secondary_chr();
		        min_primary_chr = er->get_curr_primary_chr();
                    }
                }
            }
            //}}}

	    // if the chrome pair switches, reset the max_pos
	    if ( (last_min_primary_refid != min_primary_refid) ||
		 (last_min_secondary_refid != min_secondary_refid) ) {
                max_pos = window_size;
                last_min_primary_refid = min_primary_refid;
                last_min_secondary_refid = min_secondary_refid;
            }

            bool input_processed = true;

            while (input_processed) {
                input_processed = false;
					
                //{{{ read the files, process anything in frame
                for ( i_er = inter_chrom_evidence_readers.begin();
		      i_er != inter_chrom_evidence_readers.end();
		      ++i_er) {

		    SV_EvidenceReader *er = *i_er;

		    if ( er->has_next() ) {
                        int32_t curr_primary_refid = 
                                er->get_curr_primary_refid();
			int32_t curr_secondary_refid = 
			        er->get_curr_secondary_refid();
		        CHR_POS curr_pos = er->get_curr_primary_pos();

                        if ( (curr_primary_refid <= min_primary_refid) &&
			     (curr_secondary_refid <= min_secondary_refid) &&
			     (curr_pos < max_pos) ) {

                            er->process_input_chr_pos(
                                    er->get_curr_primary_chr(),
                                    er->get_curr_secondary_chr(),
                                    max_pos,
                                    r_bin);
                            input_processed = true;
                        } 
                    }
                }
		//}}}

	        //{{{ get breakpoints
	        // get anything that has ends in both chroms
	        vector< UCSCElement<SV_BreakPoint*> > values = 
                        r_bin.values(min_secondary_chr);

                vector< UCSCElement<SV_BreakPoint*> >::iterator it;

                for (it = values.begin(); it < values.end(); ++it) {
		    SV_BreakPoint *bp = it->value;
                    if ( bp->weight >= min_weight ) {
                        //bp->do_it();
                        bp->trim_intervals();
                        bp->print_bedpe(++call_id);
                        if (show_evidence)
                            bp->print_evidence("\t");
                    }

                    if (r_bin.remove(*it, false, false, true) != 0) {
                        cerr << "Error removing element" << endl;
                        abort();
                    }
                    bp->free_evidence();
                    delete bp;
                }
                //}}}

                max_pos = max_pos * 2;
            }

	    has_next = false;
	    //{{{ Test if there is still input lines
	    for ( i_er = inter_chrom_evidence_readers.begin();
		  i_er != inter_chrom_evidence_readers.end();
		  ++i_er) {
		SV_EvidenceReader *er = *i_er;
		has_next = has_next || er->has_next();
            }
            //}}}
        }
		
        //{{{ Call remaining break points
	values = r_bin.values();

        for (it = values.begin(); it != values.end(); ++it) {
	    SV_BreakPoint *bp = it->value;
            if ( bp->weight >= min_weight ) {
	        //bp->do_it();
	        bp->trim_intervals();
                bp->print_bedpe(++call_id);
                if (show_evidence)
                    bp->print_evidence("\t");
            }

            if (r_bin.remove(*it, false, false, true) != 0) {
                cerr << "Error removing element" << endl;
                abort();
            }
            bp->free_evidence();
            delete bp;
        }
        //}}}
    }
    //}}}

    //{{{ free up stuff
    string s = inter_chrom_file_prefix + ".bam";
    remove(s.c_str());
    s = inter_chrom_file_prefix + ".sort.bam";
    remove(s.c_str());
    map<int, pair<log_space*,log_space*> >::iterator e_it;
    for( e_it =  SV_Evidence::distros.begin();
         e_it !=  SV_Evidence::distros.end();
         ++e_it) {
        free(e_it->second.first);
        free(e_it->second.second);
    }
    for ( i_er = evidence_readers.begin();
          i_er != evidence_readers.end();
          ++i_er) {
        delete(*i_er);
    }
    //}}}

    return 0;
}
