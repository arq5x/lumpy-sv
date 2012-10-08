/*****************************************************************************
 * log_space.h
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
 * ****************************************************************************/
#ifndef __LOG_SPACE_H__
#define __LOG_SPACE_H__

#include <math.h>

typedef double log_space;

double get_p(log_space ls);

log_space get_ls(double p);

// x*y
log_space ls_multiply(log_space x, log_space y);

// x/y
log_space ls_divide(log_space x, log_space y);

// x + y
log_space ls_add(log_space x, log_space y);

#endif
