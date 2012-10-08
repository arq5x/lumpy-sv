/*****************************************************************************
 * log_space.cpp
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
#include "log_space.h"
#include <iostream>

using namespace std;

double get_p(log_space ls)
{
	if ( isfinite(ls) )
		return exp(ls);
	else
		return 0;
}

log_space get_ls(double p)
{
	if (p == 0)
		return -INFINITY;
	else
		return log(p);
}

// x*y
log_space ls_multiply(log_space x, log_space y)
{
	if (isinf(x) || isinf(y))
		return -INFINITY;

	return (log_space)(x + y);
}
// x/y
log_space ls_divide(log_space x, log_space y)
{
	return (log_space)(x - y);
}

// x + y
log_space ls_add(log_space x, log_space y)
{
	if (isinf(x))
		return y;

	if (isinf(y))
		return x;

	if (x < y)
		return y + log(1 + exp(x - y));
	else
		return x + log(1 + exp(y - x));

}
