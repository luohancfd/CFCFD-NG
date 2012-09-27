// l_bc.hh

#ifndef L_BC_HH
#define L_BC_HH

#include "l_slug.hh"

int L_bc_left_velocity(GasSlug* A, double v);
int L_bc_left_reflect(GasSlug* A);
int L_bc_left_free(GasSlug* A);
int L_bc_right_velocity(GasSlug* A, double v);
int L_bc_right_reflect(GasSlug* A);
int L_bc_right_free(GasSlug* A);
int L_exchange_bc_data(GasSlug* A, GasSlug* B);
int L_blend_slug_ends(GasSlug* A, int endA, 
		      GasSlug* B, int endB,
		      double dxb);

#endif
