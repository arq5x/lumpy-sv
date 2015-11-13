/*****************************************************************************
 * SV_PairReader.cpp
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
#include "SV_Pair.h"
#include "SV_Tools.h"
#include "ucsc_bins.hpp"
#include <iostream>

//{{{ SV_PairReader:: SV_PairReader()
SV_PairReader::
SV_PairReader(struct pair_end_parameters pair_end_param)
{
    is_open = false;
    bam_file = pair_end_param.bam_file;
    sample_name = pair_end_param.sample_name;
    histo_file = pair_end_param.histo_file;
    mean = pair_end_param.mean;
    stdev = pair_end_param.stdev;
    read_length = pair_end_param.read_length;
    min_non_overlap = pair_end_param.min_non_overlap;
    discordant_z = pair_end_param.discordant_z;
    back_distance = pair_end_param.back_distance;
    weight = pair_end_param.weight;
    min_mapping_threshold = pair_end_param.min_mapping_threshold;
    ev_id = SV_EvidenceReader::counter;
    SV_EvidenceReader::counter = SV_EvidenceReader::counter + 1;
    SV_EvidenceReader::sample_names[ev_id] = sample_name;
    SV_EvidenceReader::ev_types[ev_id] = "PE";
}
//}}}

//{{{ SV_PairReader:: SV_PairReader()
SV_PairReader::
SV_PairReader()
{
    bam_file = "";
    sample_name = "";
    histo_file = "";
    mean = 0;
    stdev = 0;
    read_length = 0;
    min_non_overlap = 0;
    discordant_z = 0;
    back_distance = 0;
    weight = 0;
    min_mapping_threshold = 0;
    ev_id = SV_EvidenceReader::counter;
    SV_EvidenceReader::counter = SV_EvidenceReader::counter + 1;
    SV_EvidenceReader::sample_names[ev_id] = sample_name;
    SV_EvidenceReader::ev_types[ev_id] = "PE";

    inited = false;
}
//}}}

//{{{SV_PairReader:: ~SV_PairReader()
SV_PairReader::
~SV_PairReader()
{
    //cerr << "~SV_PairReader()" << endl;
    free(histo);
    free(distro);
}
//}}}

//{{{ string SV_PairReader:: check_params()
string
SV_PairReader::
check_params()
{
    string msg = "";

    if (bam_file.compare("") == 0)
        msg.append("bam_file ");
    if (sample_name.compare("") == 0)
        msg.append("id ");
    if (histo_file.compare("") == 0)
        msg.append("histo_file ");
    if (mean == 0)
        msg.append("mean ");
    if (stdev == 0)
        msg.append("stdev ");
    if (read_length == 0)
        msg.append("read_length ");
    if (min_non_overlap == 0)
        msg.append("min_non_overlap ");
    if (discordant_z == 0)
        msg.append("discordant_z ");
    if (back_distance == 0)
        msg.append("back_distance ");
    if (weight == 0)
        msg.append("weight ");

    return msg;
}
//}}}

//{{{ bool SV_PairReader::: add_param(string param, string val)
bool
SV_PairReader::
add_param(char *param, char *val)
{
    if ( strcmp("bam_file", param) == 0 )
        bam_file = val;
    else if ( strcmp("id", param) == 0 ) {
        sample_name = val;
	SV_EvidenceReader::sample_names[ev_id] = sample_name;
	SV_EvidenceReader::ev_types[ev_id] = "PE";
    }
    else if ( strcmp("histo_file", param) == 0 )
        histo_file = val;
    else if ( strcmp("mean", param) == 0 )
        mean = atof(val);
    else if ( strcmp("stdev", param) == 0 )
        stdev = atof(val);
    else if ( strcmp("read_length", param) == 0 )
        read_length = atoi(val);
    else if ( strcmp("min_non_overlap", param) == 0 )
        min_non_overlap = atoi(val);
    else if ( strcmp("discordant_z", param) == 0 )
        discordant_z = atoi(val);
    else if ( strcmp("back_distance", param) == 0 )
        back_distance = atoi(val);
    else if ( strcmp("weight", param) == 0 )
        weight = atoi(val);
    else if ( strcmp("min_mapping_threshold", param) == 0 )
        min_mapping_threshold = atoi(val);
    else if ( strcmp("read_group", param) == 0 )
        read_group.push_back(val);
    else
        return false;

    return true;
}
//}}}

//{{{ struct pair_end_parameters SV_PairReader:: get_pair_end_parameters()
struct pair_end_parameters
        SV_PairReader::
get_pair_end_parameters()
{
    struct pair_end_parameters pair_end_param;

    pair_end_param.bam_file = bam_file;
    pair_end_param.sample_name = sample_name;
    pair_end_param.histo_file = histo_file;
    pair_end_param.mean = mean;
    pair_end_param.stdev = stdev;
    pair_end_param.read_length = read_length;
    pair_end_param.min_non_overlap = min_non_overlap;
    pair_end_param.discordant_z = discordant_z;
    pair_end_param.back_distance = back_distance;
    pair_end_param.weight = weight;
    pair_end_param.min_mapping_threshold = min_mapping_threshold;
    pair_end_param.read_group = read_group;

    return pair_end_param;
}
//}}}

//{{{ void SV_PairReader:: set_statics()
void
SV_PairReader::
set_statics()
{
#if 0
    SV_Pair::min_mapping_threshold = min_mapping_threshold;
    SV_Pair::min_non_overlap = min_non_overlap;
    SV_Pair::insert_mean = mean;
    SV_Pair::insert_stdev = stdev;
    SV_Pair::insert_Z = discordant_z;
    SV_Pair::back_distance = back_distance;
    SV_Pair::read_length = read_length;
    SV_Pair::read_group = read_group;

    int distro_size;
    unsigned int start, end;

    SV_Pair::histo_size = read_histo_file(histo_file,
                                          &(SV_Pair::histo),
                                          &start,
                                          &end);

    SV_Pair::histo_start = start;
    SV_Pair::histo_end = end;

    SV_Pair::set_distro_from_histo();
#endif
}
//}}}

//{{{ void SV_PairReader:: unset_statics()
void
SV_PairReader::
unset_statics()
{
    //free(SV_Pair::histo);
    //free(SV_Pair::distro);
}
//}}}

//{{{ void SV_PairReader:: initialize()
void
SV_PairReader::
initialize()
{
    read_histo_file(histo_file, &histo, &histo_start, &histo_end);
    distro_size = SV_Pair::set_distro_from_histo(back_distance,
                  histo_start,
                  histo_end,
                  histo,
                  &distro);
}
//}}}

//{{{ void SV_PairReader:: process_input(
void
SV_PairReader::
process_input( BamAlignment &_bam,
               RefVector &_refs,
               UCSCBins<SV_BreakPoint*> &r_bin)
{
    if (_bam.IsMapped() && _bam.IsMateMapped())
        SV_Pair::process_pair(_bam,
                              _refs,
                              mapped_pairs,
                              r_bin,
                              weight,
                              ev_id,
                              this);
}
//}}}

//{{{ void SV_PairReader:: process_input(
void
SV_PairReader::
process_input( BamAlignment &_bam,
               RefVector &_refs,
               BamWriter &inter_chrom_reads,
               UCSCBins<SV_BreakPoint*> &r_bin)
{
    if (_bam.IsMapped() && _bam.IsMateMapped())
        SV_Pair::process_intra_chrom_pair(_bam,
                                          _refs,
                                          inter_chrom_reads,
                                          mapped_pairs,
                                          r_bin,
                                          weight,
                                          ev_id,
                                          this);
}
//}}}

//{{{ string SV_PairReader:: get_curr_chr()
string
SV_PairReader::
get_curr_chr()
{
    //cerr << "Pair get_curr_chr" << endl;
    return refs.at(bam.RefID).RefName;
}
//}}}

//{{{ CHR_POS SV_PairReader:: get_curr_pos()
CHR_POS
SV_PairReader::
get_curr_pos()
{
    //cerr << "Pair get_curr_chr" << endl;
    return bam.Position;
}
//}}}

//{{{ void SV_PairReader:: terminate()
void
SV_PairReader::
terminate()
{
    //cerr << "SplitRead Done" << endl;
    reader.Close();
}
//}}}

//{{{ bool SV_EvidenceReader:: has_next()
bool
SV_PairReader::
has_next()
{
    //cerr << "Pair has_next:" << have_next_alignment << endl;
    return have_next_alignment;
}
//}}}

//{{{ string SV_PairReader:: get_source_file_name()
string
SV_PairReader::
get_source_file_name()
{
    return bam_file;
}
//}}}

#if 0
//{{{ void SV_PairReader:: process_input_chr(string chr,
void
SV_PairReader::
process_input_chr(string chr,
                  UCSCBins<SV_BreakPoint*> &r_bin)
{
    // Process this chr, or the next chr
    while ( have_next_alignment &&
            ( chr.compare( refs.at(bam.RefID).RefName) == 0 ) ) {
        if (bam.IsMapped() && bam.IsMateMapped())  //Paired read
            SV_Pair::process_pair(bam,
                                  refs,
                                  reader.GetHeader().ToString(),
                                  mapped_pairs,
                                  r_bin,
                                  weight,
                                  ev_id,
                                  this);
        have_next_alignment = reader.GetNextAlignment(bam);
        if ( bam.RefID < 0 )
            have_next_alignment = false;
    }
}
//}}}
//{{{ void SV_PairReader:: process_input_chr(string chr,
void
SV_PairReader::
process_input_chr_pos(string chr,
                      CHR_POS pos,
                      UCSCBins<SV_BreakPoint*> &r_bin)
{
    // Process this chr, or the next chr
    while ( have_next_alignment &&
            ( chr.compare( refs.at(bam.RefID).RefName) == 0 ) &&
            ( bam.Position < pos ) ) {
        if (bam.IsMapped() && bam.IsMateMapped())  //Paired read
            SV_Pair::process_pair(bam,
                                  refs,
                                  N
                                  mapped_pairs,
                                  r_bin,
                                  weight,
                                  ev_id,
                                  this);
        have_next_alignment = reader.GetNextAlignment(bam);
        if ( bam.RefID < 0 )
            have_next_alignment = false;
    }
}
//}}}
//{{{ void SV_PairReader:: process_input()
void
SV_PairReader::
process_input( UCSCBins<SV_BreakPoint*> &r_bin)
{
    while (reader.GetNextAlignment(bam))
        if (bam.IsMapped() && bam.IsMateMapped())  //Paired read
            SV_Pair::process_pair(bam,
                                  refs,
                                  mapped_pairs,
                                  r_bin,
                                  weight,
                                  ev_id,
                                  this);
}
//}}}
#endif


