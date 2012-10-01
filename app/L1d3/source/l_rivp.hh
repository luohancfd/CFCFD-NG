// l_rivp.hh

#ifndef L_RIVP_HH
#define L_RIVP_HH

#include "l_cell.hh"

int L_rivp(const std::vector<LFlowState>& QL, 
	   const std::vector<LFlowState>& QR,
           double ustar[], double pstar[],
           int first, int last, int end1, int end2);

#endif
