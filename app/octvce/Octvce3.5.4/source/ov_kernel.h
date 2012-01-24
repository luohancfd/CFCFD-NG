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

/**\file All data structures and macros used for Octree Virtual Cell Embedding simulations*/

#include "../../Geomio/source/gio_kernel.h"
#include "../../Geomio/source/gio_adts.h"
#include "../../Geomio/source/gio_io.h"
#include "../../Geomio/source/gio_lists.h"

typedef struct flow_info Flow_info;
typedef struct state_vector State_vector;
typedef struct prim_vector Prim_vector;
typedef struct gas_dat Gas_dat;
typedef struct history_loc History_loc;
typedef struct jwlb_coeff Jwlb_coeff;
typedef struct cfl_dat Cfl_dat;
typedef struct grad Grad;
typedef struct limit Limit;
typedef struct flux_vector * Flux_vector;
typedef struct sum_flx Sum_flx;
typedef struct ic_heat IC_heat;
typedef struct deton_pts Deton_pts;
typedef struct bc_data BC_data;
typedef struct cart_cell * Cart_cell; 
typedef struct list_leaf * List_leaf; 
typedef struct list_merge * List_merge;
typedef struct list_vtx_glob * List_vtx_glob;
typedef struct point_tnode * Point_tnode;
typedef struct merge * Merge;
typedef struct vertex * Vertex;

#define NOTINCLUDE 0
#define WARHEADBURN 0
#define AXI_OK 1 /*Use 'proper' axisymmetric interface values for proper axisymmetric modelling (otherwise they're only good for no-obstruction-in-domain case)*/
#define SOLNSPLIT 1 /*Capability to resume simulation from previous solution rather than from beginning*/
#define VTXM 1 /*For vertex debugging*/
#define VTX_PAR 0 /*Parallel verticies*/
#define FLOAT_EQ_STABLE 1 /*Try a more stable form of floating point comparisons*/
#define DETON 1 /*Want to allow for finite detonations?*/
#define AFTERBURN 0
#define EXHAUSTIVE_DEBUG 0
#define CHECK_CRAP_QUANT 0
#define SUPERSONIC_VORTEX 0

/*\b SOME UBIQUITOUS MACROS*/
/*------------------------------------------------------------------*/

#define FLOAT_TOL 3e-16 /*Tolerance for floating point comparisons (slightly larger than eps)*/
#define AREA_DIFF_TOL 1e-12 /*If difference in cell areas smaller than this, don't bother setting solid boundary conditions*/
#define TRUE 1
#define FALSE 0
#define NOERROR 0
#define ERROR -1
#define ERROR2 -2
#define MAX(A,B) (((A) > (B)) ? (A) : (B))
#define MAX3(A,B,C) ((MAX(A,B) > (C)) ? (MAX(A,B)) : (C))
#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MIN3(A,B,C) ((MIN(A,B) < (C)) ? (MIN(A,B)) : (C))
#define ABS(A) (((A) < 0.0) ? (-(A)): (A))
#define ROUND(A) ((ABS((double)(A) - floor((double)(A))) >= ABS(ceil((double)(A)) - (double)(A))) ? ((int)ceil((double)(A))) : ((int)floor((double)(A))))
#define SGN(A) (((A) == 0) ? 0 : ((((double) (A))/(ABS((double) (A))))))
#define SQR(A) ((A)*(A))
#define CUBE(A) (((A)*(A))*(A))
#define QUART(A) ((((A)*(A))*(A))*(A))
#define SQRT(A) (pow(((double) (A)), 0.5))
#define EQ(A,B) (((A) >= ((B)-(FLOAT_TOL))) && ((A) <= ((B)+(FLOAT_TOL)))) /*Comparisons for stable floating point equality*/
#define NOT_EQ(A,B) (((A) < ((B)-(FLOAT_TOL))) || ((A) > ((B)+(FLOAT_TOL))))
#define UNUSED_VARIABLE(x) ((void) (x))
/*------------------------------------------------------------------*/

/*\b USEFUL MACROS FOR CELLS*/
/*------------------------------------------------------------------*/

#define MINRL 0; /**< Default minimum level of refinement of root cell*/
#define MAXRL 9; /**< Default maximum level of refinement of root cell*/

#define MAX_NUM_CHILDREN 8 /**< Obviously have 8 legitimate children for every cell in octree*/
#define MAX_NUM_NEIGHBOURS 4 /**< Each cell face must have a maximum of 4 face neighbours*/
#define NUM_FACES 6 /**< No. faces of each cell*/

#define SOLID 0 /**< Cell types for each cell*/
#define FLUID 1
#define INTERSECTED 2
#define UNKNOWN -2 /**< Many positive structure variables initialized to this*/

#define OUTSIDE -1 /**< Helps determine if a cell is inside/outside an IC volume*/
#define IRRELEVANT -3

#define ALL 0 /**< Some definitions useful for parallel adaptation*/
#define ALLOC_CHILD 3
#define POSITIVE 1
#define NEGATIVE -1
#define POSITIVE2 2
#define NEGATIVE2 -2

#define EST 0 /**< Numeric representation of 6 directions*/
#define WST 1
#define NTH 2
#define STH 3
#define UPR 4
#define LWR 5

#define CALC_CELL_EDGE_LENGTH(C) (C -> cell_length)
#define CALC_CELL_FLUID_VOL(C) (CUBE(CALC_CELL_EDGE_LENGTH(C)))
#define CALC_CELL_FLUID_AREA(C) (SQR(CALC_CELL_EDGE_LENGTH(C))) /**< Calculating basis edge length/area/volumes of cartesian cells*/

#define SMALL_CELL_LIMIT 0.05 /**< Fraction of basis cartesian volume a small cell has before it's labelled a 'small cell'*/

#define NOT_CONNECTED 0 /**< Macros for cell merging*/
#define CONNECTED 1

#define REFINE 0 /**< Adaptation status of cells after adaptation parameter(s) computed*/
#define COARSEN 1
#define REALLY_COARSEN 2 /*A child might be flagged for coarsening but be unable*/
#define JUST_MADE 4

/*------------------------------------------------------------------*/

/*Some definitions for boundary conditions*/
#define WALL 0
#define BORDER 1

#define NONREFLECTING 1
#define SPECIFIC 2
#define TRANSIENT 3
#define EXTRAPOLATED 4
#define FIXED 5

/*------------------------------------------------------------------*/

/**\b Useful equation of state definitions and macros*/

#define PERFECT_GAS 1
#define JWLB 2
#define JWLB_RHO_CUTOFF 1e-10 /*If explosion products density in cell < this value, won't apply JWLB EOS*/
#define AIR_THRESHOLD 0.004 /*Value that indicates when explosion products can be simply treated as air (if desired)*/

#define JWLB_K(I,M) (1.0 + ((M)/2.0)*(JWLB_coeffs.JRL[I][0]) +  (SQR(M)/6.0)*(JWLB_coeffs.JRL[I][1]) + (CUBE(M)/24.0)*(JWLB_coeffs.JRL[I][2]) + (QUART(M)/120.0)*(JWLB_coeffs.JRL[I][3])) /**< Related to JWLB 'q' term*/
#define JWLB_DKDM(I,M) (0.5*(JWLB_coeffs.JRL[I][0]) + ((M)/3.0)*(JWLB_coeffs.JRL[I][1]) + (SQR(M)/8.0)*(JWLB_coeffs.JRL[I][2]) + (CUBE(M)/30.0)*(JWLB_coeffs.JRL[I][3])) /**< Related to JWLB 'q' term*/

/*------------------------------------------------------------------*/

/**\b Different flux calculators*/

#define EFM 1
#define AUSM 2
#define AUSMDV 3
#define ADAPTIVE 4

/*------------------------------------------------------------------*/

/**\b STATE VECTOR*/

/*------------------------------------------------------------------*/

struct state_vector /**< \b Standard state vector for Euler euations (SI units)*/
{
  double Rhof; /*Density of gas*/
  double RhoU; 
  double RhoV; 
  double RhoW; 
  double RhoE; 
  double Pres; /*Redundant, but useful extra pressure variable*/
  double T;

  /*For the moment, only 2 species allowed (explosion products and ambient gas)*/
  double rhof_prod; /**< Density of explosion products gas*/

#if AFTERBURN
  double rhof_prod_ub; /**< Unburnt product density (for afterburning model)*/
#endif

#if SUPERSONIC_VORTEX
  double r; /*Radial distance from origin*/
#endif
};

/*------------------------------------------------------------------*/

/**\b VECTOR OF PRIMITIVE VARIABLES*/

struct prim_vector
{
  double R;
  double P;
  double T;
  double U;
  double V;
  double W; 
};

/*------------------------------------------------------------------*/

/**\b PRIMITIVE FLUX VECTOR*/

/*------------------------------------------------------------------*/

struct flux_vector /**< \b Standard generic flux vector for Euler euations (SI units).  Quantities are understood
		      to be fluxes going out of face*/
{
  double Mss_flx;
  double X_mntm_flx; 
  double Y_mntm_flx;
  double Z_mntm_flx;
  double E_flx;

  double mss_flx_prod;
#if AFTERBURN
  double mss_flx_prod_ub;
#endif
};

/*------------------------------------------------------------------*/

/**\b FLUX SUMS - STORED AREA-FLUX SUMS IF TIMESTEP NEEDS TO BE REVISED*/

/*------------------------------------------------------------------*/

struct sum_flx
{
  double sum_mss;
  double sum_xmnt;
  double sum_ymnt;
  double sum_zmnt;
  double sum_e;

  double sum_mss_prod;
#if AFTERBURN
  double sum_mss_prod_ub;
#endif
};

/*------------------------------------------------------------------*/

/**\b GAS DATA*/

struct gas_dat
{
  int num_species;

  double Cv_amb;
  double gam_amb;
    
  short int eos_prod; /*Equation of state of explosion products*/
  double Cv_prod;
  char jwlb_cutoff; /*Either 'y' for treating low densities as just air, or 'n' when JWLB EOS is always applied*/  
  double gam_prod; /*If JWLB EOS used this is irrelevant*/
};

/*------------------------------------------------------------------*/

/**\b HISTORY LOCATIONS*/

struct history_loc
{
  int num_locs;
  int num_timesteps_history;
  double **locs;
};

/*------------------------------------------------------------------*/

/**\b JWLB coefficients*/

struct jwlb_coeff
{
  short int smallest_R_index;
  short int smallest_RL_index;

  double A[5];
  double R[5];
  double AL[5];
  double BL[5];
  double RL[5];  
  double C;
  double omega;
  double v0;

  double JRL[5][4]; /**< Used in 'q' evaluation for temperature-dependent JWLB form.  RL terms by row*/
};

/*------------------------------------------------------------------*/

/**\b CFL data*/

struct cfl_dat
{
  short int CFL_ok; /**\param Helper variable to indicate if CFL no. should be revised*/
  double CFL_max; /**\param Maximum allowable CFL number*/
  double CFL_cutback; /**\param Factor by which CFL no. should be cutback*/
  double CFL; /**\param Actual minimum CFL no. (over all cells)*/
  double *CFL_min; /**\param Minimum CFL no. for each thread*/
};

/*------------------------------------------------------------------*/

/**\b HEAT ADDITION FOR EXPLOSION PRODUCTS*/

struct ic_heat
{
  int num_pts;
  int step;
  double f1, f2, e1, e2, T1, T2;
  double **heat; /**< 1st column time, 2nd is fuel consumption rate (kg/s), 3rd is energy release rate (J/(kg.s)), 
		    4th ignition temperature (K)*/
};

/*------------------------------------------------------------------*/

/**\b INFORMATION ON DETONATION*/

struct deton_pts
{
  short int deton_active; /**< Specifies if any undetonated cells still exist*/
  int num_deton_pts;
  double **dpts; /**< Each point where detonation initiated*/
  double *dcjvels; /**CJ detonation speed corresponding to each detonation point*/
};

/*------------------------------------------------------------------*/

/**\b BOUNDARY CONDITION DATA*/

struct bc_data 
{
  short int type;
  double rhof_stat[2]; /**< 0th entry indicates fixed/extrapolated, 1st entry is fixed value*/
  double pres_stat[2];
  double T_stat[2];
  double U_stat[2];
  double V_stat[2];
  double W_stat[2];

  /*If transient boundary conditions are desired*/
  short int pres_spec, rhof_spec, temp_spec;
  int num_pts;
  int step;
  char ***stat; /*Matrix containing transient values*/
  Prim_vector Pv1; /*Flow variables corresponding to transient values at RK step 1*/
  Prim_vector Pv2;
};

/*------------------------------------------------------------------*/

/**\b VECTOR OF GRADIENTS OF PRIMITIVE VARIABLES*/

struct grad
{  
  double grad_rhof[3];
  double grad_E[3];
  double grad_u[3];
  double grad_v[3];
  double grad_w[3];
};

/*------------------------------------------------------------------*/

struct limit
{
  double lim; /*Apparently 1 limiter syncs phase errors - superior convergence?*/

  double rhof_lim;
  double rhof_prod_lim;
  double E_lim;
  double u_lim;
  double v_lim;
  double w_lim;  
};

/*------------------------------------------------------------------*/

/**\b ALL REQUIRED FLOW INFO*/

/*------------------------------------------------------------------*/

struct flow_info /**< All flow variables relevant to solving flows*/
{
  State_vector State_vec;
  State_vector Temp_state_vec;

  Flux_vector Face_fluxes[NUM_FACES][MAX_NUM_NEIGHBOURS]; 
  /**< Each cell has 6 - 24 face fluxes depending on neighbours which can change.*/

  Flux_vector Obstructed_flux;
  /**< For partially obstructed cells, have pressure fluxes to account for body there*/

  Sum_flx Flx_sums; /**< If CFL needs revising, will store flux-area sum products instead of recomputation*/
  
  Grad Grads; /**< Vector of gradients for primitive variables*/
  Limit Lim; /**< Limiter values for solution reconstruction*/
};

/*------------------------------------------------------------------*/

/*\b OCTREE CELL DATA STRUCTURE*/

/*------------------------------------------------------------------*/

struct cart_cell
{
  Flow_info Flow_data; /**< Each cell has all relevant flow data*/ 

  double centroid[3]; /**< Centroid of cell in global co-ord sys*/

  Vertex verticies[MAX_NUM_CHILDREN]; /**< 8 verticies*/

  short int child_num; /**< What child (0 - 7) this cell is of its parent's; can be inferred from parent, but faster to store*/

#if SOLNSPLIT 
  char *cell_num; /*To write out grid info for later allocing and visitation, need this*/
#endif

  short int cell_level; /**< Refinement level of cell from root (level 0), can be inferred from cell no., but faster to store*/

  double cell_length; /**< Redundant, but useful edge length*/

  struct cart_cell *parent;
  
  struct cart_cell *children[MAX_NUM_CHILDREN];
  
  struct cart_cell *face_neighbours[NUM_FACES][MAX_NUM_NEIGHBOURS]; /**< At most 24 face neighbours*/
  
  double flux_area[NUM_FACES][MAX_NUM_NEIGHBOURS]; /**< Flux areas through each face quadrant*/

  double face_r_centroid[2][2]; /**< For 2D axisymmetric mode, need to know fluid interface r values (only for quadrants 0 and 1)*/

  double r_centroid; /**< For 2D axisymmetric mode, need to know cell area r value*/

  double wall_area; /**< Area of net average obstruction within a cell (can be derived from flux areas)*/

  double wall_norm[3]; /**< Normal to obstruction (wall) of cell (can be derived from flux areas) - this is into the wall*/
  
  double cell_volume; /**< Actual unobstructed volume of cell*/

  short int cell_type; /**< Cell type - FLUID, SOLID, INTERSECTED, UNKNOWN.  Is redundant given cell_volume, but
			  good for heuristics*/

  double IC_flag; /**< Specifies if cell is INSIDE/OUTSIDE/INTERSECTED by special IC region(s)*/
  
  short int un_det; /**< In finite detonation problems, a cell is un_det if a detonation wave hasn't crossed it yet*/

  List_leaf Leaf_list_loc; /**< Location on list of leaf nodes (if not on list, is NULL)*/

  Body Solid; /**< If Solid != NULL, is uniquely intersected (volume-wise) by one body*/ 
  
  short int additional_refine_times; /**< No. times cell needs additional refinement for geometrical reasons e.g. is intersecting
					a body and is too coarse etc.  If additional_refine_times > 0, cell can't be coarsened*/

  short int is_small_cell; /**< If TRUE, cell is 'small' i.e. too-obstructed; must use cell averaging*/

  Merge Merge_data; /**< All cells in merged cell would point to this*/

  List_merge Merge_list_loc; /**< Location on list of linked cells forming merged cell (stored in Merge_data)*/

  short int adjacent_to_normal_cell; /**< Identifies small cell as adjacent to normal cell or not*/

  char shocked_cell; /**< Identifies if shock is at cell - turn off reconstruction to minimize odd-even decoupling*/

  short int adapt_flag; /**< Variable to indicate adaptation status*/

  short int cell_marked; /**< Helper variable to prevent revisitation of cell during cell merging/unmerging*/

  short int just_created; /**< Helper variable indicating if this is a new cell on the list*/

  List_bbox Intersected_bbox; /**< Temporary list of intersected bodies (candidates)*/

  short int bbox_unknown; /**< Helper flag for Intersected_bbox*/
};

/*------------------------------------------------------------------*/

/**\b LIST DATA STRUCTURES*/

/*------------------------------------------------------------------*/

struct list_leaf /**< \b Leaf nodes of octree - purely to avoid expensive octree traversals*/
{
  Cart_cell cell_loc; /**< Address of leaf cell*/
  short int pure_fluid; /**< Helper variable during mesh and flow soln writing - identifies FLUID cells next to wall*/
  int thread_num; /**< Denotes which thread is processing this list node*/
  
  struct list_leaf *prev; struct list_leaf *next; 
};

struct list_merge /**< \b List of cells linked to each other to form a merged cell*/
{
  Cart_cell cell_loc; 

  struct list_merge *prev; struct list_merge *next;
};

struct list_vtx_glob /**< \b Node on global list of verticies*/
{
  Vertex vtx_loc;
  short int delete; /*Mark a vertex for deletion from the list*/

  struct list_vtx_glob *prev; struct list_vtx_glob *next;
};

/*------------------------------------------------------------------*/

/**\b POINT ADT DATA STRUCTURE*/

/*------------------------------------------------------------------*/

struct point_tnode /**< \b ADT of points in 2D/3D*/
{
  double point[3]; /**< If 2D space, point[2] redundant*/
  int level; /**< Level of node on ADT*/

  double lower_bound[3];
  double upper_bound[3];
  /**< Partitioned region corresponding to point in 2D/3D space (if 2D space, 3rd co-ordinate redundant)*/

  struct point_tnode *left; struct point_tnode *right;
  
  int point_num; /**< point_num th body on ADT*/
  struct point_tnode *parent;
};

/*------------------------------------------------------------------*/

/**\b DATA ON MERGED CELLS*/

/*------------------------------------------------------------------*/

struct merge
{
  double cell_volume;
  List_merge Linked_cells; /*List of cells linked together to form merged cells*/
};

/*------------------------------------------------------------------*/

/**\b VERTEX DATA STRUCTURE*/

/*------------------------------------------------------------------*/

struct vertex
{
  double loc[3]; /*Location of vertex*/
  int *vtx_nums;  /*Vertex number fo plotting mesh (could be different for different threads)*/
  
  Cart_cell Leaf_cells[MAX_NUM_CHILDREN]; /*Stores leaf cells sharing this vertex (cell position in array depends on
					    its spatial relationship to vertex) - max 8 cells can share vertex*/

  List_vtx_glob Vtxlist_loc; /*Position on global list of verticies*/
};

/*------------------------------------------------------------------*/
