/*****************************************************************************
 * SV_InterChromBamReader.cpp
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
#include "SV_InterChromBamReader.h"
#include "SV_Pair.h"
#include "SV_SplitRead.h"
#include "SV_Tools.h"
#include "ucsc_bins.hpp"
#include <iostream>
#include <map>

//{{{ SV_InterChromBamReader:: SV_InterChromBamReader()
SV_InterChromBamReader::
SV_InterChromBamReader()
{
}
//}}}

//{{{ SV_InterChromBamReader:: SV_InterChromBamReader(string _inter_chrom_bam..
SV_InterChromBamReader::
SV_InterChromBamReader(string _inter_chrom_bam_file,
		       map<pair<string,string>, SV_EvidenceReader*> *_bam_evidence_readers)
{
	inter_chrom_bam_file = _inter_chrom_bam_file,
	bam_evidence_readers = _bam_evidence_readers;
}
//}}}

//{{{SV_InterChromBamReader:: ~SV_InterChromBamReader()
SV_InterChromBamReader::
~SV_InterChromBamReader()
{
    /*
    map<pair<string,string>, SV_EvidenceReader*>::iterator it;
    for (it = bam_evidence_readers->begin();
            it != bam_evidence_readers->end();
            ++it)
        delete(it->second);
    */
}
//}}}

//{{{ string SV_InterChromBamReader:: check_params()
string
SV_InterChromBamReader::
check_params()
{
	return "";
}
//}}}

//{{{ bool SV_InterChromBamReader::: add_param(string param, string val)
bool
SV_InterChromBamReader::
add_param(char *param, char *val)
{
	return false;
}
//}}}

//{{{ void SV_InterChromBamReader:: initialize()
void
SV_InterChromBamReader::
initialize()
{

	if (!bam_reader.Open(inter_chrom_bam_file)) {
		cerr << "Could not open temp file" << inter_chrom_bam_file << 
				" for writing... Aborting." << endl;
		abort();
	}

	refs = bam_reader.GetReferenceData();
	header = bam_reader.GetHeader().ToString();
	has_next_alignment = bam_reader.GetNextAlignment(bam);
        bam.QueryBases.clear();
        bam.AlignedBases.clear();
        bam.Qualities.clear();
}
//}}}

//{{{ bool SV_InterChromBamReader:: has_next()
bool
SV_InterChromBamReader::
has_next()
{
	return has_next_alignment;
}
//}}}

#if 1
//{{{ string SV_InterChromBamReader:: get_curr_chr()
string
SV_InterChromBamReader::
get_curr_chr()
{
	//cerr << "Pair get_curr_chr" << endl;
	return refs.at(bam.RefID).RefName;
}
//}}}
#endif

//{{{ string SV_InterChromBamReader:: get_curr_primary_refid()
int32_t
SV_InterChromBamReader::
get_curr_primary_refid()
{
	//if (bam.IsFirstMate())
	if (bam.RefID < bam.MateRefID)
		return bam.RefID;
	else
		return bam.MateRefID;
}
//}}}

//{{{ string SV_InterChromBamReader:: get_curr_secondary_refid()
int32_t
SV_InterChromBamReader::
get_curr_secondary_refid()
{
	if (bam.RefID < bam.MateRefID)
		return bam.MateRefID;
	else
		return bam.RefID;
}
//}}}

//{{{ string SV_InterChromBamReader:: get_curr_primary_chr()
string
SV_InterChromBamReader::
get_curr_primary_chr()
{
	//if (bam.IsFirstMate())
	if (bam.RefID < bam.MateRefID)
		return refs.at(bam.RefID).RefName;
	else
		return refs.at(bam.MateRefID).RefName;
}
//}}}

//{{{ string SV_InterChromBamReader:: get_curr_secondary_chr()
string
SV_InterChromBamReader::
get_curr_secondary_chr()
{
	//if (bam.IsFirstMate())
	if (bam.RefID < bam.MateRefID)
		return refs.at(bam.MateRefID).RefName;
	else
		return refs.at(bam.RefID).RefName;
}
//}}}

//{{{ CHR_POS SV_InterChromBamReader:: get_curr_primary_pos()
CHR_POS
SV_InterChromBamReader::
get_curr_primary_pos()
{
	if (bam.IsFirstMate())
		return bam.Position;
	else
		return bam.MatePosition;

}
//}}}

//{{{ CHR_POS SV_InterChromBamReader:: get_curr_primary_pos()
CHR_POS
SV_InterChromBamReader::
get_curr_secondary_pos()
{
	if (bam.IsFirstMate())
		return bam.MatePosition;
	else
		return bam.Position;

}
//}}}

//{{{ void SV_InterChromBamReader:: set_statics()
void
SV_InterChromBamReader::
set_statics()
{
	// the statics will be set while the input is being proccessed
}
//}}}

#if 0
//{{{ void SV_InterChromBamReader:: process_input()
void
SV_InterChromBamReader::
process_input( UCSCBins<SV_BreakPoint*> &r_bin)
{
	exit(1);
}
//}}}
#endif

#if 0
//{{{ void SV_InterChromBamReader:: process_input(
void
SV_InterChromBamReader::
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
//{{{ void SV_InterChromBamReader:: process_input_chr(string chr,
void
SV_InterChromBamReader::
process_input_chr(string chr,
				  UCSCBins<SV_BreakPoint*> &r_bin)
{
	exit(1);
}
//}}}
#endif

//{{{ void SV_InterChromBamReader:: process_input_chr_pos(string chr,
void
SV_InterChromBamReader::
process_input_chr_pos(string primary_chr,
					  string secondary_chr,
					  CHR_POS pos,
					  UCSCBins<SV_BreakPoint*> &r_bin)
{
	// Process this chr, or the next chr 
	while ( has_next_alignment &&
			( primary_chr.compare( get_curr_primary_chr() ) == 0 ) &&
			( secondary_chr.compare( get_curr_secondary_chr() ) == 0 ) &&
			( get_curr_primary_pos() < pos ) ) {

		//get the name of the current reader from the LS tag
		string file_name;
		bam.GetTag("LS",file_name);

		// get the name of the current read_group from the RG tag
		string curr_read_group = "";
		bam.GetTag("RG",curr_read_group);
		pair<string,string> ev_pair (file_name,curr_read_group);

		SV_EvidenceReader *curr_reader;
		// first search for matching bamfile,read_group SV_EvidenceReader
		map<pair<string,string>, SV_EvidenceReader*>::iterator it;
		it = (*bam_evidence_readers).find(ev_pair);
		if (it != (*bam_evidence_readers).end()) {
		    curr_reader = it->second;
		    curr_reader->process_input(bam,refs,r_bin);
		}
		else {
		  // then search for matching bamfile SV_EvidenceReader
		  it = (*bam_evidence_readers).find(pair<string,string> (file_name,""));
		  if (it != (*bam_evidence_readers).end()) {
		      curr_reader = it->second;
		      curr_reader->process_input(bam,refs,r_bin);
		  }
		}

		has_next_alignment = bam_reader.GetNextAlignment(bam);
                bam.QueryBases.clear();
                bam.AlignedBases.clear();
                bam.Qualities.clear();
#if 0
		curr_reader = (*bam_evidence_readers)[bam.Filename];
		last_file = bam.Filename;
		curr_reader->process_input(bam,refs,inter_chrom_reads,r_bin);
		last_reader = curr_reader;

		has_next_alignment = bam_reader.GetNextAlignment(bam);

		if ( bam.RefID < 0 )
			has_next_alignment = false;
#endif
	}
}
//}}}

//{{{ void SV_InterChromBamReader:: terminate()
void 
SV_InterChromBamReader::
terminate()
{
	bam_reader.Close();
}
//}}}

//{{{ string SV_InterChromBamReader:: get_source_file_name()
string
SV_InterChromBamReader::
get_source_file_name()
{
	return inter_chrom_bam_file;
}
//}}}
