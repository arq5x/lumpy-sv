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
class SV_BreakPoint;
#ifndef __SV_Evidence_H__
#define __SV_Evidence_H__

#include <iostream>
#include <map>
#include <utility>

#include "SV_BreakPoint.h"
#include "SV_Evidence.h"
#include "SV_EvidenceReader.h"

using namespace std;

#include "ucsc_bins.hpp"
#include "log_space.h"
#include <map>

class SV_Evidence
{

public:
        static map<int, pair<log_space*,log_space*> > distros;
        static map<int, int> distros_size;
        static UCSCBins<int> exclude_regions;
        int weight;
        int ev_id;
        int type;

        virtual void set_bp_interval_probability(struct breakpoint_interval *i);

        virtual string evidence_type();

        virtual void print_evidence();
        virtual ~SV_Evidence();
        virtual SV_BreakPoint* get_bp();

};
#endif
