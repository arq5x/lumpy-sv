/*****************************************************************************
 * SV_EvidenceReader.cpp
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

#include "SV_EvidenceReader.h"
#include "api/BamReader.h"
using namespace BamTools;

using namespace std;

int  SV_EvidenceReader:: counter = 1;
map<int,string> SV_EvidenceReader:: sample_names;
map<int,string> SV_EvidenceReader:: ev_types;

//{{{ SV_EvidenceReader:: ~SV_EvidenceReader()
SV_EvidenceReader::
~SV_EvidenceReader()
{
    //cerr << "~SV_EvidenceReader()" << endl;
}
//}}}

//{{{ string SV_EvidenceReader:: check_params()
string
SV_EvidenceReader::
check_params()
{
	string msg = "";

	return msg;
}
//}}}

//{{{ bool SV_EvidenceReader::: add_param(string param, string val)
bool
SV_EvidenceReader::
add_param(char *param, char *val)
{
	return false;
}
//}}}

//{{{ void SV_EvidenceReader:: initialize()
void
SV_EvidenceReader::
initialize()
{
}
//}}}

//{{{ void SV_EvidenceReader:: set_statics()
void
SV_EvidenceReader::
set_statics()
{
}
//}}}

//{{{ void SV_EvidenceReader:: unset_statics()
void
SV_EvidenceReader::
unset_statics()
{
}
//}}}

#if 0
//{{{ void SV_EvidenceReader:: process_input( BamAlignment &_bam,
void
SV_EvidenceReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
			   string header,
			   UCSCBins<SV_BreakPoint*> &r_bin)
{
}
//}}}

//{{{ void SV_EvidenceReader:: process_input(UCSCBins<SV_BreakPoint*> &l_bin,
void
SV_EvidenceReader::
process_input( UCSCBins<SV_BreakPoint*> &r_bin)
{
}
//}}}

//{{{ void SV_EvidenceReader:: process_input_chr(string chr,
void
SV_EvidenceReader::
process_input_chr(string chr,
				  UCSCBins<SV_BreakPoint*> &r_bin)
{
}
//}}}
#endif

//{{{ void SV_EvidenceReader:: process_input_chr_pos(string chr,
void
SV_EvidenceReader::
process_input_chr_pos(string chr,
					  CHR_POS pos,
					  UCSCBins<SV_BreakPoint*> &r_bin)
{
}
//}}}

//{{{void SV_EvidenceReader:: process_input( BamAlignment &_bam,
void
SV_EvidenceReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
			   BamWriter &inter_chrom_reads,
			   UCSCBins<SV_BreakPoint*> &r_bin)
{
	abort();
}
//}}}

//{{{void SV_EvidenceReader:: process_input( BamAlignment &_bam,
void
SV_EvidenceReader::
process_input( BamAlignment &_bam,
			   RefVector &_refs,
			   UCSCBins<SV_BreakPoint*> &r_bin)
{
	abort();
}
//}}}


//{{{ void SV_EvidenceReader:: terminate()
void
SV_EvidenceReader::
terminate()
{
}
//}}}

//{{{ string SV_EvidenceReader:: get_curr_chr()
string
SV_EvidenceReader::
get_curr_chr()
{
	return "";
}
//}}}

//{{{ string SV_EvidenceReader:: get_curr_pos()
CHR_POS
SV_EvidenceReader::
get_curr_pos()
{
	CHR_POS x = 0;
	return x;
}
//}}}

//{{{ bool SV_EvidenceReader:: has_next()
bool
SV_EvidenceReader::
has_next()
{
	return false;
}
//}}}

//{{{ string SV_EvidenceReader:: get_source_file_name()
string
SV_EvidenceReader::
get_source_file_name()
{
	return "Error";
}
//}}}
//

//{{{ string SV_EvidenceReader:: get_curr_primary_refid()
int32_t
SV_EvidenceReader::
get_curr_primary_refid()
{
	cerr << "Error reaching SV_EvidenceReader:: get_curr_primary_refid "
		<<endl;
	abort();
	int32_t a = 1;
	return a;
}
//}}}

//{{{ string SV_EvidenceReader:: get_curr_secondary_refid()
int32_t
SV_EvidenceReader::
get_curr_secondary_refid()
{
	cerr << "Error reaching SV_EvidenceReader:: get_curr_secondary_refid "
		<<endl;
	abort();
	int32_t a = 1;
	return a;
}
//}}}


//{{{ string SV_EvidenceReader:: get_curr_primary_chr()
string
SV_EvidenceReader::
get_curr_primary_chr()
{
	cerr << "Error reaching SV_EvidenceReader:: get_curr_primary_chr "
		<<endl;
	abort();
	string a = "";
	return a;
}
//}}}

//{{{ string SV_EvidenceReader:: get_curr_secondary_chr()
string
SV_EvidenceReader::
get_curr_secondary_chr()
{
	cerr << "Error reaching SV_EvidenceReader:: get_curr_secondary_chr "
		<<endl;
	abort();
	string a = "";
	return a;
}
//}}}

//{{{ CHR_POS SV_EvidenceReader:: get_curr_primary_pos()
CHR_POS
SV_EvidenceReader::
get_curr_primary_pos()
{
	abort();
	CHR_POS a = 0;
	return a;
}
//}}}

//{{{ CHR_POS SV_EvidenceReader:: get_curr_primary_pos()
CHR_POS
SV_EvidenceReader::
get_curr_secondary_pos()
{
	abort();
	CHR_POS a = 0;
	return a;
}
//}}}

//{{{ void SV_EvidenceReader:: process_input_chr_pos(string chr,
void
SV_EvidenceReader::
process_input_chr_pos(string primary_chr,
					  string secondary_chr,
					  CHR_POS pos,
					  UCSCBins<SV_BreakPoint*> &r_bin)
{
	cerr << "Error reaching SV_EvidenceReader:: process_input_chr_pos "
		<<endl;
	abort();
}
//}}}
