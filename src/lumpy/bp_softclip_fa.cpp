/*****************************************************************************
  bamToBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "version.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BamAncillary.h"
#include "bedFile.h"
#include "sequenceUtils.h"
#include "SV_Evidence.h"
#include "SV_Pair.h"
#include "SV_BreakPoint.h"

#include "ucsc_bins.hpp"

using namespace BamTools;

#include <vector>
#include <map>
#include <algorithm>    // for swap()
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <gsl_statistics_int.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bamToBedpe"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


//{{{ forward declarations
void ShowHelp(void);
//}}}

//{{{ void ShowHelp(void) {
void ShowHelp(void) {
    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" <<
			endl <<
		"Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl <<
		"Summary: Converts BAM alignments to BED6 or BEDPE format." << endl <<
			endl <<
		"Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl <<
		"Options: " << endl <<
		"\t-novo\t" << "Pairs were aligned using Novoalign" << endl <<
					   endl <<
		"\t-bwa\t"  << "Pairs were aligned using BWA" << endl <<
					   endl <<
		 "\t-x\t"   << "The minimum amount of non-overlap between" << endl <<
					   "\t\tthe two mapped reads in a pair." << endl <<
					   "\t\tDefault is 10." << endl <<
					   endl <<
		"\t-m\t"	<< "Insert mean" << endl <<
					   endl <<
		"\t-s\t"	<< "Insert stdev" << endl <<
					   endl <<
		"\t-z\t"	<< "Insert Z" << endl <<
					   endl <<
		"\t-tt\t"	<< "Trim threshold" << endl <<
					   endl <<
		"\t-mt\t"	<< "Merge threshold" << endl <<
					   endl <<
		"\t-me\t"	<< "Minimum evidence" << endl <<
					   endl <<
		"\t-cz\t"	<< "Cluser Z (Default 2)" << endl <<
					   endl <<
		"\t-mno\t"	<< "Minimum non-overlap(Default 20)" << endl <<
					   endl <<
		"\t-sz\t"	<< "Evidence size Z (Default 10)" << endl <<
					   endl <<
		endl;
	
    // end the program here
    exit(1);
}
//}}}

//{{{ int main(int argc, char* argv[]) {
int main(int argc, char* argv[])
{

	//{{{ setup
    // our configuration variables
    string bamFile = "stdin";
    string tag = "";

    bool haveBam  = false;
    bool useNovoalign = false;
    bool useBWA = false;
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
    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
				haveBam = true;
                bamFile = argv[i + 1];
                i++;
            }
        }

        else if(PARAMETER_CHECK("-novo", 5, parameterLength))
                useNovoalign = true;

        else if(PARAMETER_CHECK("-bwa", 4, parameterLength)) 
                useBWA = true;

		/*
        else if(PARAMETER_CHECK("-z", 2, parameterLength)) {
            if ((i+1) < argc) {
				haveZ = true;
                Z = atoi(argv[i + 1]);
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
	//}}}

    //{{{ make sure we have an input files and either novo or bwa
    if (haveBam == false) {
        cerr << endl << "*****" << endl << 
				"*****ERROR: Need -i (BAM) file. " << endl << "*****" << endl;
        ShowHelp();
    }
    if (useNovoalign == false && useBWA == false) {
	    cerr << endl << "*****" << endl <<
				"*****ERROR: Must select either -novo or -bwa" << 
				endl << "*****" << endl;
        ShowHelp();

	}

	//}}}
	
	SV_Pair::min_non_overlap = min_non_overlap;
	SV_Pair::insert_mean = mean;
	SV_Pair::insert_stdev = stdev;
	SV_Pair::insert_Z = Z;
	SV_Pair::cluster_Z = cluster_Z;

	SV_BreakPoint::p_trim_threshold = trim_threshold;
	SV_BreakPoint::p_merge_threshold = merge_threshold;

	UCSCBins<SV_BreakPoint*> bins;

    // open the BAM file
    BamReader reader;
    reader.Open(bamFile);

    // get header & reference information
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

	// make map for matching pairs
	map<string, BamAlignment> mapped_pairs;
	map<string, BamAlignment> orphan_pairs;

    // rip through the BAM file and convert each mapped entry to BED
    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {
		if (bam.IsMapped() && bam.IsMateMapped()) { //Paired read
			mapped_pair(bam, refs, mapped_pairs, bins);
		}
		// split reads
		// orphan pairs
	}
    reader.Close();


	vector< UCSCElement<SV_BreakPoint*> > values = bins.values();
	vector< UCSCElement<SV_BreakPoint*> >::iterator it;

	sort(values.begin(), values.end(),
			UCSCElement<SV_BreakPoint*>::sort_ucscelement_by_value);

	it = unique(values.begin(),
				values.end(),
				UCSCElement<SV_BreakPoint*>::compare_ucscelement_by_value);

	values.resize( it - values.begin() );

	sort(values.begin(), values.end(),
			UCSCElement<SV_BreakPoint*>::sort_ucscelement_by_start);

	unsigned int id = 0;

	int *sizes = (int *) malloc (values.size() * sizeof(int));
	int i = 0;

	for (it = values.begin(); it < values.end(); ++it) {
		SV_BreakPoint *bp = it->value;
		if ( bp->evidence.size() > min_evidence) {
			sizes[i] = bp->evidence.size();
			++i;
		}
	}

	double bp_evidence_count_mean = gsl_stats_int_mean(sizes, 1, i);
	double bp_evidence_count_stdev = gsl_stats_int_sd_m(sizes, 1, i, 
			bp_evidence_count_mean);

	cerr << bp_evidence_count_mean << "\t" << bp_evidence_count_stdev << endl;

	for (it = values.begin(); it < values.end(); ++it) {
		SV_BreakPoint *bp = it->value;
		if ( (bp->evidence.size() >= min_evidence) &&
			 (bp->evidence.size() <= (
					bp_evidence_count_mean + 
					bp_evidence_count_stdev*bp_evidence_count_Z) ) ) {
			 
			bp->trim_intervals();
			bp->print_bedpe();
			//cout << *bp << endl;

			/*
			vector<SV_Evidence*>::iterator it;
			for (it = bp->evidence.begin(); it < bp->evidence.end(); ++it) {
				SV_Evidence *sv_e = *it;
				cout << "\t";
				sv_e->print_evidence();
				cout << endl;
			}
			*/
		}
	}

}
//}}}
