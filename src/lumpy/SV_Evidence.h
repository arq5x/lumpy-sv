/*****************************************************************************
 * SV_Evidence.h
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

#ifndef __SV_Evidence_H__
#define __SV_Evidence_H__

#include <iostream>

using namespace std;

#include "log_space.h"

class SV_Evidence
{

	public:
		virtual void print_evidence();
		int weight;
		int id;
		int type;
		virtual ~SV_Evidence();
};
#endif
