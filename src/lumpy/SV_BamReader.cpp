/*****************************************************************************
 * SV_BamReader.cpp
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

#include "SV_PairReader.h"
#include "SV_SplitReadReader.h"
#include "SV_BamReader.h"
#include "SV_Pair.h"
#include "SV_SplitRead.h"
#include "SV_Tools.h"
#include "ucsc_bins.hpp"
#include <iostream>
#include <map>

//{{{ SV_BamReader:: SV_BamReader()
SV_BamReader::
SV_BamReader()
{
	curr_reader = NULL;
	tmp_file_name = "/tmp/lumpy_inter_chrom_reads.bam";
}
//}}}

//{{{ SV_BamReader:: SV_BamReader(map<string, SV_EvidenceReader*>
SV_BamReader::
SV_BamReader(map<string, SV_EvidenceReader*> *_bam_evidence_readers)
{
	bam_evidence_readers = _bam_evidence_readers;
	curr_reader = NULL;
	tmp_file_name = "/tmp/lumpy_inter_chrom_reads.bam";
}
//}}}

//{{{SV_BamReader:: ~SV_BamReader()
SV_BamReader::
~SV_BamReader()
{
	map<string, SV_EvidenceReader*>::iterator it;
	for (it = bam_evidence_readers->begin();
			it != bam_evidence_readers->end();
			++it)
		delete(it->second);
}
//}}}

//{{{ string SV_BamReader:: check_params()
string
SV_BamReader::
check_params()
{
	return "";
}
//}}}

//{{{ bool SV_BamReader::: add_param(string param, string val)
bool
SV_BamReader::
add_param(char *param, char *val)
{
	return false;
}
//}}}

//{{{ void SV_BamReader:: initialize()
void
SV_BamReader::
initialize()
{
	vector<string> bam_files;
	map<string, SV_EvidenceReader*>::iterator it;
	for (it = bam_evidence_readers->begin(); 
		 it != bam_evidence_readers->end(); 
		 ++it) 
		bam_files.push_back(it->first);

	if ( !bam_reader.Open(bam_files) ) {
		cerr << "Could not open input BAM files." << endl;
		exit(1);
	}    
	refs = bam_reader.GetReferenceData();
	header = bam_reader.GetHeader().ToString();
	has_next_alignment = bam_reader.GetNextAlignment(bam);

	if ( !inter_chrom_reads.Open(tmp_file_name,
								 header,
								 refs) ) {
		cerr << "Could not open temp file" << tmp_file_name << 
				" for writing... Aborting." << endl;
		abort();
	}
}
//}}}

//{{{ bool SV_BamReader:: has_next()
bool
SV_BamReader::
has_next()
{
	return has_next_alignment;
}
//}}}

//{{{ string SV_BamReader:: get_curr_chr()
string
SV_BamReader::
get_curr_chr()
{
	//cerr << "Pair get_curr_chr" << endl;
	return refs.at(bam.RefID).RefName;
}
//}}}

//{{{ CHR_POS SV_BamReader:: get_curr_pos()
CHR_POS
SV_BamReader::
get_curr_pos()
{
	return bam.Position;
}
//}}}

//{{{ void SV_BamReader:: set_statics()
void
SV_BamReader::
set_statics()
{
	// the statics will be set while the input is being proccessed
}
//}}}

#if 0
//{{{ void SV_BamReader:: process_input()
void
SV_BamReader::
process_input( UCSCBins<SV_BreakPoint*> &r_bin)
{
	exit(1);
}
//}}}
#endif

#if 0
//{{{ void SV_BamReader:: process_input(
void
SV_BamReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
			   string header,
			   UCSCBins<SV_BreakPoint*> &r_bin)
{
	exit(1);
}
//}}}
#endif

#if 0
//{{{ void SV_BamReader:: process_input_chr(string chr,
void
SV_BamReader::
process_input_chr(string chr,
				  UCSCBins<SV_BreakPoint*> &r_bin)
{
	exit(1);
}
//}}}
#endif

//{{{ void SV_BamReader:: process_input_chr_pos(string chr,
void
SV_BamReader::
process_input_chr_pos(string chr,
					  CHR_POS pos,
					  UCSCBins<SV_BreakPoint*> &r_bin)
{
	string last_file = "";
	SV_EvidenceReader *last_reader;

	// Process this chr, or the next chr 
	while ( has_next_alignment &&
			( chr.compare( refs.at(bam.RefID).RefName) == 0 ) &&
			( bam.Position < pos ) ) {
		curr_reader = (*bam_evidence_readers)[bam.Filename];
		last_file = bam.Filename;
		curr_reader->process_input(bam,refs,inter_chrom_reads,r_bin);
		last_reader = curr_reader;

		has_next_alignment = bam_reader.GetNextAlignment(bam);

		if ( bam.RefID < 0 )
			has_next_alignment = false;
	}
}
//}}}

//{{{ void SV_BamReader:: terminate()
void 
SV_BamReader::
terminate()
{
	bam_reader.Close();
	inter_chrom_reads.Close();
}
//}}}

//{{{ string SV_BamReader:: get_source_file_name()
string
SV_BamReader::
get_source_file_name()
{
	return "";
}
//}}}

//{{{ void SV_BamReader:: set_inter_chrom_file_name()
void 
SV_BamReader::
set_inter_chrom_file_name(string file_name)
{
	tmp_file_name = file_name;
}
//}}}
