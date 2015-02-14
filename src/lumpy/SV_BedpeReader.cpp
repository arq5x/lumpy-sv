/*****************************************************************************
 * SV_BedpeReader.cpp
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

#include "SV_Bedpe.h"
#include "SV_BedpeReader.h"
#include "SV_Tools.h"
#include <stdlib.h>
#include <iostream>

//{{{ SV_BedpeReader:: SV_BedpeReader()
SV_BedpeReader::
SV_BedpeReader()
{
    bedpe_file = "";
    sample_name = "";
    //distro_file = "";
    distro_size = 0;
    weight = 0;
    ev_id = SV_EvidenceReader::counter;
    SV_EvidenceReader::counter = SV_EvidenceReader::counter + 1;
    SV_EvidenceReader::sample_names[ev_id] = sample_name;
    SV_EvidenceReader::ev_types[ev_id] = "BD";
}
//}}}

//{{{ string SV_BedpeReader:: check_params()
string
SV_BedpeReader::
check_params()
{
    string msg = "";

    if (bedpe_file.compare("") == 0)
        msg.append("bedpe_file ");
    if (sample_name.compare("") == 0)
        msg.append("id ");
    /*
    if (distro_file.compare("") == 0)
        msg.append("distro_file ");
    */

    if (weight == 0)
        msg.append("weight ");

    return msg;
}
//}}}

//{{{ bool SV_BedpeReader::: add_param(string param, string val)
bool
SV_BedpeReader::
add_param(char *param, char *val)
{


    if ( strcmp("bedpe_file", param) == 0 )
        bedpe_file = val;
    else if ( strcmp("id", param) == 0 ) {
        sample_name = val;
	SV_EvidenceReader::sample_names[ev_id] = sample_name;
	SV_EvidenceReader::ev_types[ev_id] = "BD";
    }
    /*
    else if ( strcmp("distro_file", param) == 0 )
        distro_file = val;
    */
    else if ( strcmp("weight", param) == 0 )
        weight = atoi(val);
    else
        return false;

    return true;
}
//}}}

//{{{ struct pair_end_parameters SV_BedpeReader:: get_pair_end_parameters()
struct bedpe_parameters
        SV_BedpeReader::
get_bedpe_parameters()
{
    struct bedpe_parameters bedpe_param;

    bedpe_param.bedpe_file = bedpe_file;
    bedpe_param.sample_name = sample_name;
    //bedpe_param.distro_file = distro_file;
    bedpe_param.weight = weight;

    return bedpe_param;
}
//}}}

//{{{ void SV_BedpeReader:: set_statics()
void
SV_BedpeReader::
set_statics()
{
#if 0
    SV_Bedpe::distro_size = read_distro_file(distro_file,
                            &(SV_Bedpe::distro),
                            &(SV_Bedpe::distro_start),
                            &(SV_Bedpe::distro_end));

    SV_Bedpe::back_distance = back_distance;
#endif
}
//}}}

//{{{ void SV_BedpeReader:: unset_statics()
void
SV_BedpeReader::
unset_statics()
{
    //free(SV_Bedpe::distro);
}
//}}}

//{{{ void SV_BedpeReader:: initialize()
void
SV_BedpeReader::
initialize()
{
    /*
    distro_size = read_distro_file(distro_file,
                                   &distro,
                                   &distro_start,
                                   &distro_end);
    */
    bedpe= new BedFilePE(bedpe_file);

    lineNum = 0;
    bedpe->Open();
    bedpeStatus = bedpe->GetNextBedPE(bedpeEntry, lineNum);

}
//}}}

//{{{ void SV_BedpeReader:: process_input()
void
SV_BedpeReader::
process_input( UCSCBins<SV_BreakPoint*> &r_bin)
{
    while (	bedpeStatus != BED_INVALID) {
        if (bedpeStatus == BED_VALID) {
            SV_Bedpe::process_bedpe(&bedpeEntry,
                                    r_bin,
                                    weight,
                                    ev_id,
                                    this);
            bedpeEntry = nullBedpe;
        }

        bedpeStatus = bedpe->GetNextBedPE(bedpeEntry, lineNum);
    }
}
//}}}

//{{{ string SV_BedpeReader:: get_curr_chr()
string
SV_BedpeReader::
get_curr_chr()
{
    return bedpeEntry.chrom1;
}
//}}}

//{{{ string SV_BedpeReader:: get_curr_chr()
CHR_POS
SV_BedpeReader::
get_curr_pos()
{
    if (bedpeEntry.start1 < bedpeEntry.start2)
        return bedpeEntry.start1;
    else
        return bedpeEntry.start2;
}
//}}}

//{{{ void SV_BedpeReader:: process_input_chr(string chr,
void
SV_BedpeReader::
process_input_chr(string chr,
                  UCSCBins<SV_BreakPoint*> &r_bin)
{

    while ( ( bedpeStatus != BED_INVALID ) &&
            ( chr.compare( get_curr_chr() ) == 0 ) ) {


        if (bedpeStatus == BED_VALID) {
            SV_Bedpe::process_bedpe(&bedpeEntry,
                                    r_bin,
                                    weight,
                                    ev_id,
                                    this);
            bedpeEntry = nullBedpe;
        }

        bedpeStatus = bedpe->GetNextBedPE(bedpeEntry, lineNum);
    }
}
//}}}

//{{{ void SV_BedpeReader:: process_input_chr_pos(string chr,
void
SV_BedpeReader::
process_input_chr_pos(string chr,
                      CHR_POS pos,
                      UCSCBins<SV_BreakPoint*> &r_bin)
{

    while ( ( bedpeStatus != BED_INVALID ) &&
            ( chr.compare( get_curr_chr() ) == 0 ) &&
            ( bedpeEntry.start1  < pos ) &&
            ( bedpeEntry.start2  < pos ) ) {

        if (bedpeStatus == BED_VALID) {
            SV_Bedpe::process_bedpe(&bedpeEntry,
                                    r_bin,
                                    weight,
                                    ev_id,
                                    this);
            bedpeEntry = nullBedpe;
        }

        bedpeStatus = bedpe->GetNextBedPE(bedpeEntry, lineNum);
    }
}
//}}}

//{{{ void SV_BedpeReader:: terminate()
void
SV_BedpeReader::
terminate()
{
    bedpe->Close();
}
//}}}

//{{{ bool SV_EvidenceReader:: has_next()
bool
SV_BedpeReader::
has_next()
{
    return  (bedpeStatus != BED_INVALID);
}
//}}}
