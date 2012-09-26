// l_misc.hh

#ifndef L_MISC_HH
#define L_MISC_HH

#include "l_kernel.hh"
#include "l_tube.hh"

int L_alloc(struct slug_data *A);
void L_free(struct slug_data *A);
int L_set_index_range(struct slug_data *A);
int L_copy_cell_data(struct L_cell *source,
                     struct L_cell *target, int copy_extras);
int L_interpolate_cell_data(struct slug_data *A, 
			    double xloc, struct L_cell& icell);
double L_slug_end_pressure(struct slug_data *A, int which_end, double dx);
int L_slug_end_properties(struct slug_data *A, int which_end, double dx, 
			  double * total_mass, struct L_flow_state * Q);
int L_blend_cells( struct L_cell *cA, struct L_cell *cB, 
		   struct L_cell *c, double alpha, int blend_type );
int L_blend_slug_ends(struct slug_data *A, int endA, 
		      struct slug_data *B, int endB,
		      double dxb );
int L_fill_data(struct slug_data *A);
int L_compute_areas(struct slug_data *A, TubeModel *tube);
int L_dump_cell(struct slug_data *A, int ix);
int maximum_p(struct slug_data *A, double *p_max, double *x_max);
int total_energy(struct slug_data *A, double *E_tot);
double L_get_dt_plot(SimulationData *SD);
double L_get_dt_history(SimulationData *SD);

#endif
