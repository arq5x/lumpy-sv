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

using namespace std;

//{{{ SV_EvidenceReader:: ~SV_EvidenceReader()
//SV_EvidenceReader::
//~SV_EvidenceReader()
//{
//}
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

//{{{ void SV_EvidenceReader:: process_input(UCSCBins<SV_BreakPoint*> &l_bin,
void
SV_EvidenceReader::
process_input(UCSCBins<SV_BreakPoint*> &l_bin,
						   UCSCBins<SV_BreakPoint*> &r_bin)
{
}
//}}}

//{{{ void SV_EvidenceReader:: process_input_chr(string chr,
void
SV_EvidenceReader::
process_input_chr(string chr,
				  UCSCBins<SV_BreakPoint*> &l_bin,
				  UCSCBins<SV_BreakPoint*> &r_bin)
{
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

//{{{ bool SV_EvidenceReader:: has_next()
bool
SV_EvidenceReader::
has_next()
{
	return false;
}
//}}}
