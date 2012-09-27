// l_rivp.hh

#ifndef L_RIVP_HH
#define L_RIVP_HH

#include "l_tstep.hh"

int L_rivp(struct L_flow_state QL[], struct L_flow_state QR[],
           double ustar[], double pstar[],
           int first, int last, int end1, int end2);

#endif
