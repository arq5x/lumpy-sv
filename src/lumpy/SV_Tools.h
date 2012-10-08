/*****************************************************************************
 * SV_Tools.h
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

#ifndef __SV_TOOLS_H__
#define __SV_TOOLS_H__

#include <string>
using namespace std;

int read_histo_file(string file_name,
					 double **histo,
					 unsigned int *start,
					 unsigned int *end);

int read_distro_file(string file_name,
					 double **distro,
				     int *start,
				     int *end);

#endif
