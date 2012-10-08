/*****************************************************************************
 * SV_SplitReadReader.cpp
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

#include "SV_SplitReadReader.h"
#include "SV_SplitRead.h"
#include <stdlib.h>
#include <iostream>

//{{{ SV_SplitReadReader:: SV_SplitReadReader()
SV_SplitReadReader::
SV_SplitReadReader()
{
	bam_file = "";
	back_distance = 0;
	min_non_overlap = 25;
	weight = -1;
	min_mapping_threshold = 0;
	id = -1;
}
//}}}

//{{{ string SV_SplitReadReader:: check_params()
string
SV_SplitReadReader::
check_params()
{
	string msg = "";

	if (bam_file.compare("") == 0)
		msg.append("bam_file ");
	if (back_distance == 0)
		msg.append("back_distance ");
	if (weight == 0)
		msg.append("weight ");
	if (id == -1)
		msg.append("id ");

	return msg;
}
//}}}

//{{{ bool SV_SplitReadReader::: add_param(string param, string val)
bool
SV_SplitReadReader::
add_param(char *param, char *val)
{
	if ( strcmp("bam_file", param) == 0 )
		bam_file = val;
	else if ( strcmp("min_non_overlap", param) == 0 )
		min_non_overlap = atoi(val);
	else if ( strcmp("back_distance", param) == 0 )
		back_distance = atoi(val);
	else if ( strcmp("weight", param) == 0 )
		weight = atoi(val);
	else if ( strcmp("id", param) == 0 )
		id = atoi(val);
	else if ( strcmp("min_mapping_threshold", param) == 0 )
		min_mapping_threshold = atoi(val);
	else 
		return false;

	return true;
}
//}}}

//{{{ void SV_SplitReadReader:: set_statics()
void
SV_SplitReadReader::
set_statics()
{
	SV_SplitRead::min_mapping_threshold = min_mapping_threshold;
	SV_SplitRead::back_distance = back_distance;
	SV_SplitRead::min_non_overlap = min_non_overlap;
	//SV_SplitRead:: read_length = read_length;
	//SV_SplitRead:: min_split_size = min_split_size;
}
//}}}

//{{{ void SV_SplitReadReader:: initialize()
void
SV_SplitReadReader::
initialize()
{
	//cerr << "SplitRead initialize" << endl;
	// open the BAM file
	reader.Open(bam_file);

	// get header & reference information
	header = reader.GetHeaderText();
	refs = reader.GetReferenceData();

	have_next_alignment = reader.GetNextAlignment(bam);
	//cerr << "SplitRead have_next_alignment " << have_next_alignment << endl;
}
//}}}

//{{{ string SV_SplitReadReader:: get_curr_chr()
string
SV_SplitReadReader::
get_curr_chr()
{
	return refs.at(bam.RefID).RefName;
}
//}}}

//{{{ void SV_SplitReadReader:: process_input()
void
SV_SplitReadReader::
process_input(UCSCBins<SV_BreakPoint*> &l_bin,
			  UCSCBins<SV_BreakPoint*> &r_bin)
{
	while (reader.GetNextAlignment(bam)) 
		SV_SplitRead::process_split(bam,
									refs,
									mapped_splits,
									l_bin,
									r_bin,
									weight,
									id);
}
//}}}

//{{{ void SV_SplitReadReader:: process_input_chr(string chr,
void
SV_SplitReadReader::
process_input_chr(string chr,
				  UCSCBins<SV_BreakPoint*> &l_bin,
				  UCSCBins<SV_BreakPoint*> &r_bin)
{
	cerr << "SplitRead:" << chr << endl;
	// Process this chr, or the next chr 
	while ( have_next_alignment &&
			( chr.compare( refs.at(bam.RefID).RefName) == 0 ) ) {

		SV_SplitRead::process_split(bam,
									refs,
									mapped_splits,
									l_bin,
									r_bin,
									weight,
									id);

		have_next_alignment = reader.GetNextAlignment(bam);
		if ( bam.RefID < 0 )
			have_next_alignment = false;
	}
}
//}}}

//{{{ void SV_SplitReadReader:: terminate()
void 
SV_SplitReadReader::
terminate()
{
	//cerr << "SplitRead Done" << endl;
	reader.Close();
}
//}}}

//{{{ bool SV_EvidenceReader:: has_next()
bool
SV_SplitReadReader::
has_next()
{
	return have_next_alignment;
}
//}}}
