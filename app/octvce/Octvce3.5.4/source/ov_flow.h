/*Copyright (C) Joseph Tang, 2007
                                                                                                        
    This file is part of OctVCE.
                                                                                                        
    OctVCE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
                                                                                                        
    OctVCE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
                                                                                                        
    You should have received a copy of the GNU General Public License
    along with OctVCE.  If not, see <http://www.gnu.org/licenses/>.
*/

/**\file Prototypes to solve Euler equations and flow related quantities*/

short int march_soln(void);

int count_cells(void);

void get_global_timestep(void);

short int get_flow_soln(double *);

short int time_integrate(List_leaf, List_leaf, short int, int, short int);

short int all_cells_exchange_fluxes(List_leaf, List_leaf, double, short int);

void compute_fluxes(State_vector *, State_vector *, double [], Flux_vector, short int);

double min_local_cell_timestep(Cart_cell);

int exxef(double, double *, double *);

#if SUPERSONIC_VORTEX
double sv_dens(double, double, double, double);
short int calc_vortex_err(double, double, double, int, int, double);
#endif
