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

/**\file This is where all the global variables are kept*/

#include <stdio.h>
#include <stdlib.h>
#include "ov_kernel.h"

int stepnow;
int RK; /**< Runge-Kutta step*/

/**\param 2D parameters*/
short int twoD = 0; /**< Default is 3D code - if '1' is 2D, '2' is 2D axisymmetric*/
double twoD_thickness; /**< Gap thickness (which should be half length of highest refined cell)*/

/**\param Root cell*/
Cart_cell Root = NULL;

/**\param Globally accessible ADTs (pointers to root) for traversal (for cells that need to interrogate geometry engine)*/
Point_tnode ADT2Dpoint = NULL; /**< ADT of points in 2D (to determine partial areas)*/
Point_tnode ADT3Dpoint = NULL; /**< ADT of points in 3D (to determine partial volumes)*/

/**\param List of leaf nodes of octree (NULL to begin with) and list of all verticies*/
List_leaf Leaves = NULL;
List_leaf Leaves_tail = NULL;
List_vtx_glob Vtxlist = NULL;
List_vtx_glob Vtxlist_tail = NULL;
List_leaf *Leaf_heads; /**< Will contain an array of pointers to heads of leaf cell sublists*/  
List_leaf *Leaf_tails; /**< Pointers to list tails*/
List_vtx_glob *Vtx_heads;
List_vtx_glob *Vtx_tails;
double * min_dt; /**< Array of minimum timesteps*/
int *num_cells_per_thread; /**< No. of cells for each thread*/

short int output_grid = FALSE; /**< If we want to output the whole grid and solution to continue it at a later time*/
short int output_grid_last = FALSE; /**< If we only want to output the whole grid at the very last timestep*/

/**\param Global variables that control degree and nature of adaptation*/
short int num_times_refine_root; /**< No. times root is to be refined during initial grid generation*/
short int min_refinement_level; /**< Minimum level of refinement of root cell*/
short int max_refinement_level; /**< Maximum level of refinement of root cell*/
short int min_intersect_refine_level; /**< Minimum level of refinement an intersected cell*/
short int IC_level; /**< Minimum level of refinement of cell that intersects IC volume*/

short int want_to_adapt = TRUE; /**< Turn adaptation on/off*/
short int adapt_type = 1; /**< '1' uses velocity differences, '2' uses density derivatives*/
double adapt_param = 0.05; /**< Adaptation threshold for compression tolerance*/
double adapt_noise = 0.05; /**< Noise filter value for indicator type 2*/
double refine_param = 0.3; /**< Threshold for refinement (Petrie's criterion)*/
double coarsen_param = 0.1; /**< Threshold for coarsening (Petrie's criterion)*/

int num_timesteps_adapt = 10; /**< Adapt every set no. of time steps*/

/**\param Global variables that deal with spatial properties of root cell*/
double root_shift[3] = {0, 0, 0}; /**< Translation factors (shift of root's centroid from origin of global co-ord system)*/
double root_scale = 1.0; /**< Scale factor of root cell's edge lengths compared to unity*/

/**\param Type of flux solver*/
short int flux_type = UNKNOWN;

/**\param On limiting and solution order*/
char many_limit;
char always_many_limit;
double switch_order = 0.0;

/**\param Gas properties*/
Gas_dat gas_data;

/**\param Maximum density gradient (for schlieren normalization)*/
double max_schl;

/**\param History files*/
History_loc Histories;
char base_hist_filename[128];
  
/**\param JWLB parameters for explosion products gas*/
Jwlb_coeff JWLB_coeffs;

/**\param CFL data*/
Cfl_dat CFL_data;

/**\param Parameters related to timestepping*/
short int set_global_dt = FALSE; /*Fixes global timestep*/
double global_dt; /*Global minimum timestep*/ 
int num_timesteps_run = 0; /*Run simulation for how many timesteps?*/
double finish_time = 0;
double current_flow_time = 0;

/**\param No. threads for parallel processing*/
int number_of_threads = 1;
short int adapt_in_parallel = FALSE;

/**\param Variables to specify how geometry engine is to be used*/
short int use_geom_engine_always  = TRUE; /**< Interrogate geometry engine at each stage of cell refinement?  TRUE/FALSE*/

short int level_before_can_interrogate = 0; /**< Interrogate geometry engine (use VCE) after root refined 
					       level_before_can_interrogate_times.  Ensures too-coarse cells won't return erroneous 
					       results if VCE used on them.*/ 

short int refining_root_now = UNKNOWN; /**< Let's geometry engine know if root is currently being refined (initially)*/

short int allow_multiple_refinements = FALSE; /**< Allow geometry engine to tell refine() that a cell should be refined more than
						 once because of intersection reasons (TRUE/FALSE)*/

short int area_subcell_num = 20; /**< No. subcells along a face edge to determine partial areas (must be even)*/

short int volume_subcell_num = 10; /**< No. subcells along a cell edge to determine partial volumes (should be even)*/

short int wall_rep = 1; /**< How to treat wall fluxes - either in 'staircased' or 'smoothed' fashion (using the surface normal)*/

/**\param State vectors of initial conditions*/
State_vector IC_ambient;
State_vector **IC_states = NULL; /**< An array describing all IC states*/
Body *IC_regions = NULL; /**< Array of IC regions*/
double ***IC_bboxes; /**< Array of bounding boxes in 3D*/
int num_ICs = 0; /**< No. of specific IC conditions (except for ambient)*/
Deton_pts Deton; /**< Stores information on detonation locations*/

/**\param Transient heat addition for initial explosion products*/
IC_heat ICheat;

/**\param Boundary data on border faces of root cell*/
BC_data BC_est;
BC_data BC_wst;
BC_data BC_nth;
BC_data BC_sth;
BC_data BC_upr;
BC_data BC_lwr;

/**\param Some visualization parameters*/
short int only_get_boundary_cells_flow = FALSE; /**< When writing flow soln, only display boundary cells or not*/
int num_timesteps_dump_data = 100; /**< How many timesteps before writing flow solns, meshes*/
double time_interval_write = 0;
char base_soln_filename[128];

#if SUPERSONIC_VORTEX
double M_in = 4.5;
double max_drr_tol = 1e-10;
#endif
