/*****************************************************************************
 * SV_Tools.h
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
 * ***************************************************************************/

#ifndef __SV_TOOLS_H__
#define __SV_TOOLS_H__

#include <string>
#include "ucsc_bins.hpp"

#include "log_space.h"
using namespace std;

//static inline int strnum_cmp(const char *a, const char *b);

int read_histo_file(string file_name,
		    double **histo,
                    unsigned int *start,
                    unsigned int *end);

int read_distro_file(string file_name,
	             double **distro,
                     int *start,
                     int *end);

bool sort_inter_chrom_bam(string in_file_name,
			  string out_file_name);

bool create_sorted_temp_file(vector<BamAlignment>& buffer,
                             string out_file_name,
                             int num_runs,
                             string header_text,
                             RefVector &ref);

bool merge_sorted_files(string out_file_name,
			int buff_count,
                        string header_text,
                        RefVector &ref);

bool write_temp_file(vector<BamAlignment>& buffer,
		     string temp_file_name,
                     string header_text,
                     RefVector &ref);


void normalize_ls(CHR_POS size, log_space *o, log_space *r);

void parse_exclude_file(string exclude_bed_file,
                        UCSCBins<int> &exclude_regions);
uint32_t
count_clipped(vector< CigarOp > cigar_data);

#endif
