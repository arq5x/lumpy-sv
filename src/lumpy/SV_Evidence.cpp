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
	for (j = 0; j < size; ++j) {
		tmp_p[j] = src_p[j];
	}

	i->p = tmp_p;
}
//}}}
