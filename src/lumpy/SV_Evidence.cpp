/*****************************************************************************
 * SV_Evidence.cpp
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

#include "SV_Evidence.h"
#include "SV_BreakPoint.h"
#include <iostream>

using namespace std;

map<int, pair<log_space*,log_space*> > SV_Evidence:: distros;
map<int, int> SV_Evidence:: distros_size;
UCSCBins<int> SV_Evidence:: exclude_regions;

void
SV_Evidence::
print_evidence()
{
    cerr << "***** :( *****" << endl;
}

SV_Evidence::
~SV_Evidence()
{
    //cerr << "~SV_Evidence() ***** :( *****" << endl;
}

SV_BreakPoint*
SV_Evidence::
get_bp()
{
    return NULL;
}

string
SV_Evidence::
evidence_type()
{
    return "";
}
//{{{ void SV_Pair:: set_interval_probability()
void
SV_Evidence::
set_bp_interval_probability(struct breakpoint_interval *i)
{
    int size = i->i.end - i->i.start + 1;
    log_space *tmp_p = (log_space *) malloc(size * sizeof(log_space));
    log_space *src_p;

    unsigned int j;
    if (i->i.strand == '+')
        src_p = SV_Evidence::distros[ev_id].first;
    else
        src_p = SV_Evidence::distros[ev_id].second;


    // It is possible that the start of this interval was trucated because it
    // started close the the start of the chrome, and the back distance for the
    // + strand or the extension of the distribution for the - strand would
    // have caused an underrun.  In this case we need to clip the begining of
    // the distro
    // It is also possible that the end of this inteval was trucated because it
    // overlapped its pair's read interval
    // the full size is from 0 to size
    // i->i.start_clip and i->i.end_clip tell us how much of each side of the
    // disto to remove

    // if src_p is NULL then we assume a uniform distribution

    int offset = i->i.start_clip;
    for (j = 0; j < size; ++j) {
        tmp_p[j] = src_p[j + offset];
    }

    i->p = tmp_p;
}
//}}}
