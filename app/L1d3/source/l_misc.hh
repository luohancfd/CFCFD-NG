// l_misc.hh

#ifndef L_MISC_HH
#define L_MISC_HH

#include "l_kernel.hh"
#include "l_tube.hh"

int L_interpolate_cell_data(struct slug_data *A, 
			    double xloc, struct L_cell& icell);
double L_slug_end_pressure(struct slug_data *A, int which_end, double dx);
int L_slug_end_properties(struct slug_data *A, int which_end, double dx, 
			  double * total_mass, struct L_flow_state * Q);
int L_blend_slug_ends(struct slug_data *A, int endA, 
		      struct slug_data *B, int endB,
		      double dxb );

#endif
