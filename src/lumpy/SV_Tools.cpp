/*****************************************************************************
 * SV_Tools.cpp
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

#include "SV_Tools.h"

#include <string>

using namespace std;

//{{{int read_histo_file(string file_name,
/*
 * the histo file contains the histogram (frequency) of insert sizes
 * in the format
 * i	freq_i
 * i+1	freq_i+1
 * ...
 * n	freq_n
 *
 * start will be set to i 
 * and will be set to n
 *
 * i..n is assumed to be continuous (no gaps)
 */
int read_histo_file(string file_name,
					 double **histo,
					 unsigned int *start,
					 unsigned int *end)
{
	int file_len = 0;
	char line[LINE_MAX];
	FILE *file = fopen(file_name.c_str(), "r");

	if (file == NULL) {
		fprintf(stderr, "Could not open %s\n", file_name.c_str());
	}

	while ( fgets(line, LINE_MAX, file) )
		++file_len;

	fclose(file);

	*histo = (double *) malloc(file_len * sizeof(double));
	file = fopen(file_name.c_str(), "r");
	int i = 0;
	bool first = true;
	while ( fgets(line, LINE_MAX, file) ) {
		// we only need to read start, and from that value we can infer the
		// value of end
		unsigned int tmp = atoi(strtok(line, "\t"));
		if (first == true) {
			*start = tmp;
			first = false;
		}
		(*histo)[i]= atof(strtok(NULL, "\t"));
		++i;
	}

	*end = *start + file_len;

	fclose(file);

	return file_len;
}
//}}}

//{{{ int read_distro_file(string file_name,
int read_distro_file(string file_name,
					 double **distro,
				     int *start,
				     int *end)
{
	int file_len = 0;
	char line[LINE_MAX];
	FILE *file = fopen(file_name.c_str(), "r");

	if (file == NULL) {
		fprintf(stderr, "Could not open %s\n", file_name.c_str());
	}

	while ( fgets(line, LINE_MAX, file) )
		++file_len;

	fclose(file);

	*distro = (double*) malloc(file_len * sizeof(double));

	file = fopen(file_name.c_str(), "r");
	int i = 0;
	bool first = true;
	int last;
	while ( fgets(line, LINE_MAX, file) ) {
		int v = atof(strtok(line, "\t"));
		if (first == true) {
			*start = v;
			first = false;
		}
		last = v;
		(*distro)[i]= atof(strtok(NULL, "\t"));
		++i;
	}

	*end = last;

	fclose(file);

	return file_len;
}
//}}}

