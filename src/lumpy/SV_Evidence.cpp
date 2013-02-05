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

void
SV_Evidence::
print_evidence()
{
	cerr << "***** :( *****";
}

SV_Evidence::
~SV_Evidence()
{
}

SV_BreakPoint*
SV_Evidence::
get_bp()
{
	return NULL;
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
		src_p = SV_Evidence::distros[sample_id].first;
	else
		src_p = SV_Evidence::distros[sample_id].second;

	// It is possible that the start of this interval was trucated because it
	// started close the the start of the chrome, and the back distance for the
	// + strand or the extension of the distribution for the - strand would
	// have caused an underrun.  In this case we need to clip the begining of
	// the distro
	int offset = SV_Evidence::distros_size[sample_id] - size;
	for (j = 0; j < size; ++j) {
		tmp_p[j] = src_p[j + offset];
	}

	i->p = tmp_p;
}
//}}}
