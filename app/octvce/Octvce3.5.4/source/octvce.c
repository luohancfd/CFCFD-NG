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

/*Test the flow solver*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ov_kernel.h"
#include "ov_adts.h"
#include "ov_adapt.h"
#include "ov_flow.h"
#include "ov_io.h"
#include "time_keeper.h"

#define NOTINCLUDE 0
#define SPECIAL 0
#define SPECIAL2 0
#define NORMAL 1

extern State_vector **IC_states;

extern short int twoD;
extern Cart_cell Root;
extern List_leaf Leaves;
extern List_leaf Leaves_tail;
extern List_leaf * Leaf_heads;
extern List_leaf * Leaf_tails;
extern List_vtx_glob * Vtx_heads;
extern List_vtx_glob * Vtx_tails;
extern double *min_dt;
extern int * num_cells_per_thread;
extern Point_tnode ADT2Dpoint;
extern Point_tnode ADT3Dpoint;
extern short int min_intersect_refine_level;
extern short int num_times_refine_root;
extern short int min_refinement_level;
extern short int max_refinement_level;
extern short int min_intersect_refine_level;
extern short int IC_level;
extern short int refining_root_now;
extern short int set_global_dt;
extern State_vector IC_ambient;
extern int num_ICs;
extern Deton_pts Deton;
extern Gas_dat gas_data;
extern Jwlb_coeff JWLB_coeffs;
extern Cfl_dat CFL_data;
extern double root_scale;
extern int number_of_threads;
extern IC_heat ICheat;
extern short int adapt_in_parallel;
extern short int output_grid;
extern short int output_grid_last;
extern double switch_order;

extern double root_3D_bbox[2][3];
extern Body_tnode ADTbody;

int main(int argc, char *argv[])
{
  Cart_cell C;
  List_leaf L,Dummy_list; 
  int i;
  int num_bodies = 0;
  short int ic_set = FALSE;
  short int bc_set = FALSE;
  short int par_set = FALSE;
  short int gas_set = FALSE;
  short int deton_set = FALSE;
  short int geom_set = FALSE;
  short int count_energy = FALSE;
  short int continue_soln = FALSE;
  char icfile[256] = {"\0"}, bcfile[256] = {"\0"}, parfile[256] = {"\0"}, gasfile[256] = {"\0"}, detonfile[256] = {"\0"};
  char geom_file[256] = {"\0"}, grid_file[256] = {"\0"};
  double *ec, *vc;
  int *nc;
				 
  /*For some puzzling reason, I should initialize some global variables here as ISO C90 doesn't like it when I do it elsewhere*/
  num_times_refine_root = MINRL;
  min_refinement_level = MINRL;
  max_refinement_level = MAXRL;
  min_intersect_refine_level = MINRL;
  IC_level = MINRL;
  ICheat.num_pts = 0; /*Initialize this first (in case heat addition is unnecessary)*/
  Deton.num_deton_pts = 0; /*Initialize this first (in case instantaneous detonation is sufficient)*/
  Deton.deton_active = FALSE;

  gas_data.num_species = 1;  
  
#if NOTINCLUDE
  number_of_threads = 1;
  read_IC("ov_IC.par");
  read_BC("ov_BC.par"); 
  read_ov_param("ov.par");
  read_ov_gas("ov_gas.par");
  num_bodies = read_body_files("bodies.pvtp.xml");
  ic_set = TRUE; bc_set = TRUE; par_set = TRUE;
#else 
  for(i = 0; i < argc; i++) {
    if(strcmp(argv[i], "-ncpus") == 0) { /*Set no. threads/cpus for parallel solution*/
      if((i+1) >= argc) {
	printf("octvce() - No value for no. cpus.  Exiting ...\n");
	return(ERROR);
      }
      else {
	number_of_threads = (int) strtol(argv[i+1], NULL, 10);
	if(number_of_threads <= 0) {
	  printf("octvce() - No. cpus given is \"%s\", which is in error.  Exiting ...\n", argv[i+1]);
	  return(ERROR);
	} 
#ifdef _OPENMP
	if(omp_get_num_procs() == 1) {
	  printf("octvce() - Warning, only 1 physical processor seemingly available on this system.  Will set no. threads to 1.\n");
	  number_of_threads = 1; /*Do we want this step?*/
	}
#else
	printf("octvce() - Warning.  Not compiled with OpenMP enabled, so no. threads/cpus set to 1 anyway.\n");
	number_of_threads = 1;
#endif
      }      
    }
    else if(strcmp(argv[i], "-switch-order") == 0) { /*Some time you want to switch solution order?*/
      if((i+1) >= argc) {
	printf("octvce() - No time specified at which order should be switched.  Exiting ...\n");
	return(ERROR);
      }
      else {
	switch_order = strtod(argv[i+1], NULL);
	if(switch_order <= 0.0) {
	  printf("octvce() - Time for switching of order <= 0.  Don't want to switch order?  Then don't include argument.  Exiting ...\n");
	  return(ERROR);
	}
      }
    }
    else if(strcmp(argv[i], "-gas") == 0) {
      if((i+1) >= argc) {
	printf("octvce() - No file specified for gas properties.  Exiting ...\n");
	return(ERROR);
      }
      else strcpy(gasfile, argv[i+1]);
      gas_set = TRUE;
    }
    else if(strcmp(argv[i], "-ic") == 0) {
      if((i+1) >= argc) {
	printf("octvce() - No file specified for initial condition data.  Exiting ...\n");
	return(ERROR);
      }
      else strcpy(icfile, argv[i+1]); 
      ic_set = TRUE;
    }
    else if(strcmp(argv[i], "-bc") == 0) {
      if((i+1) >= argc) {
	printf("octvce() - No file specified for boundary condition data.  Exiting ...\n");
	return(ERROR);
      }
      else strcpy(bcfile, argv[i+1]);
      bc_set = TRUE;      
    }
    else if(strcmp(argv[i], "-par") == 0) {
      if((i+1) >= argc) {
	printf("octvce() - No file specified for OctVCE parameter data.  Exiting ...\n");
	return(ERROR);
      }
      else strcpy(parfile, argv[i+1]);
      par_set = TRUE;
    }
    else if(strcmp(argv[i], "-geom") == 0) {
      if((i+1) >= argc) {
	printf("octvce() - No VTK pvtp file specified for geometry data.  Exiting ...\n");
	return(ERROR);
      }
      else strcpy(geom_file, argv[i+1]);
      geom_set = TRUE;
    }
    else if(strcmp(argv[i], "-deton") == 0) {
      if((i+1) >= argc) {
	printf("octvce() - No file specified for detonation data.  Exiting ...\n");
	return(ERROR);
      }
      else strcpy(detonfile, argv[i+1]);
      deton_set = TRUE; 
    }
    else if(strcmp(argv[i], "-2D=axisymmetric") == 0) {
      printf("Now entering 2D axisymmetric mode.  Please ensure gap thickness = 0.5 of length of highest refined cell, and\n");
      printf("please ensure you are working in XY plane and X is symmetry axis and Y is radial coordinate\n");
      twoD = 2;
    }
    else if(strcmp(argv[i], "-2D") == 0) {
      printf("Now entering 2D mode.  Please ensure gap thickness = 0.5 of length of highest refined cell, and\n");
      printf("please ensure you are working in XY plane\n");
      twoD = 1;
    }
    else if(strcmp(argv[i], "-count=energy") == 0) {
      count_energy = TRUE;
    }
    else if(strcmp(argv[i], "-continue-soln") == 0) {
      if((i+1) >= argc) {
	printf("octvce() - No file specified for grid and mesh data for solution continuation.  Exiting ...\n");
	return(ERROR);
      }
      else strcpy(grid_file, argv[i+1]);
      continue_soln = TRUE;
    }
    else if(strcmp(argv[i], "-output-grid") == 0) {
      output_grid = TRUE;
    }
    else if(strcmp(argv[i], "-output-grid=last") == 0) {
      output_grid = TRUE;
      output_grid_last = TRUE;
    }
  }
#endif

  if((ic_set == FALSE) && (continue_soln == FALSE)) {
    printf("octvce() - Didn't see initial condition data file in command line arguments (\"-ic <file>\").  Exiting ...\n");
    return(ERROR);
  }  
  if(bc_set == FALSE) {
    printf("octvce() - Didn't see boundary condition data file in command line arguments (\"-bc <file>\").  Exiting ...\n");
    return(ERROR);
  }
  if(par_set == FALSE) {
    printf("octvce() - Didn't see OctVCE parameter file in command line arguments (\"-par <file>\").  Exiting ...\n");
    return(ERROR);
  }
  if(gas_set == FALSE) {
    printf("octvce() - Didn't see gas parameter file in command line arguments (\"-gas <file>\").  Exiting ...\n");
    return(ERROR);
  }

  /*Read these parameters in the order below as some later functions can be switched off (from read_ov_param())*/
  
  if(geom_set == TRUE) {
    num_bodies = read_body_files(geom_file); 
    if(num_bodies == ERROR) 
      return(ERROR);
  }

  if(read_ov_gas(gasfile) == ERROR)
    return(ERROR);

  if(read_ov_param(parfile) == ERROR)
    return(ERROR);

  if((twoD != 0) && (adapt_in_parallel == TRUE)) {
    printf("octvce(): Warning - in 2D modes, can't adapt in parallel because I can't be bothered.  But you still get parallel flow solution and I/O\n");
    adapt_in_parallel = FALSE;
  }  

  if(continue_soln == FALSE)
    if(read_IC(icfile) == ERROR)
      return(ERROR);
  
  if(deton_set == TRUE)
    if(read_deton(detonfile) == ERROR)
      return(ERROR);

  if(read_BC(bcfile) == ERROR)
    return(ERROR);
  
#if 0
  for(i = 0; i < 5; i++)
    printf("JWLB.A[%hd] = %g\n",i,JWLB_coeffs.A[i]);
  for(i = 0; i < 5; i++)
    printf("JWLB.R[%hd] = %g\n",i,JWLB_coeffs.R[i]);
  for(i = 0; i < 5; i++)
    printf("JWLB.AL[%hd] = %g\n",i,JWLB_coeffs.AL[i]);
  for(i = 0; i < 5; i++)
    printf("JWLB.BL[%hd] = %g\n",i,JWLB_coeffs.BL[i]);
  for(i = 0; i < 5; i++)
    printf("JWLB.RL[%hd] = %g\n",i,JWLB_coeffs.RL[i]);  
  printf("JWLB.C = %g\n",JWLB_coeffs.C);
  printf("JWLB.omega = %g\n",JWLB_coeffs.omega);
  printf("JWLB.v0 = %g\n",JWLB_coeffs.v0);

  printf("TNT initial soundspd is %g, squared %g, P %e\n", get_mixture_soundspd(IC_states[0]), pow(get_mixture_soundspd(IC_states[0]),2),
	 IC_states[0]->Pres);
  printf("TNT initial energy is %g\n",get_mixture_RhoE(IC_states[0])/(IC_states[0]->Rhof));
  abort();
#endif

  /*If ncpus isn't set, no. threads will be 1, if no geometry read in, will assume no geometry*/

  Root = make_root();
  if(Root == NULL) {
    printf("octvce() - Allocating root cell unsuccessful.  Exiting ...\n");
    return(ERROR);
  }

  for(i = 0; i < 3; i++)
    {
#if VTXM
      root_3D_bbox[0][i] = Root -> verticies[0] -> loc[i];
      root_3D_bbox[1][i] = Root -> verticies[7] -> loc[i];
#endif
#if !VTXM
      root_3D_bbox[0][i] = Root -> centroid[i] - 0.5*root_scale;
      root_3D_bbox[1][i] = Root -> centroid[i] + 0.5*root_scale;
#endif
    }

  if(num_bodies <= 0) 
    {
      printf("octvce(): Warning - read_body_files() couldn't see any geometries.  Will assume no geometry.\n");

      /*Root automatically has FLUID interior*/

      Root -> cell_type = FLUID;
      Root -> cell_volume = CALC_CELL_FLUID_VOL(Root);
      for(i = 0; i < NUM_FACES; i++) {
	Root -> flux_area[i][0] = CALC_CELL_FLUID_AREA(Root)/4.0;
	Root -> flux_area[i][1] = CALC_CELL_FLUID_AREA(Root)/4.0;
	Root -> flux_area[i][2] = CALC_CELL_FLUID_AREA(Root)/4.0;
	Root -> flux_area[i][3] = CALC_CELL_FLUID_AREA(Root)/4.0;
      }
      Root -> additional_refine_times = 0;
    }

  if(build_volume_ADT() == ERROR) {
    printf("octvce() - Building volume ADT unsuccessful.  Exiting ...\n");
    return(ERROR);
  }
      
  if(build_area_ADT() == ERROR) {
    printf("octvce() - Building area ADT unsuccessful.  Exiting ...\n");
    return(ERROR);
  }

  Dummy_list = NULL; 

  min_dt = malloc(number_of_threads*sizeof(double));
  (CFL_data.CFL_min) = malloc(number_of_threads*sizeof(double));
  CFL_data.CFL_cutback = 0.1; /*Initialize to 0.1 for now*/
  CFL_data.CFL_ok = TRUE; /*CFL no. is OK for now*/

  if(set_global_dt == TRUE)
    CFL_data.CFL = 1.0;
  else CFL_data.CFL = CFL_data.CFL_max; /*Is maximum CFL for now*/

  /*Set no. threads for parallel processing*/
  if(number_of_threads > 1)
    {
#ifdef _OPENMP
      omp_set_num_threads(number_of_threads);
#endif
      /*Allocate pointers to sublists (as number of threads aren't expected to change)*/
      Leaf_heads = malloc(number_of_threads*sizeof(List_leaf)); /*For heads of different sublists for each thread*/
      Leaf_tails = malloc(number_of_threads*sizeof(List_leaf)); /*Get tails of sublists*/
      Vtx_heads = malloc(number_of_threads*sizeof(List_vtx_glob));
      Vtx_tails = malloc(number_of_threads*sizeof(List_vtx_glob));

      num_cells_per_thread = malloc(number_of_threads*sizeof(int));      
    }      
  
  if(continue_soln == TRUE) {
    printf("octvce() - now reading in previous solution and grid file '%s'\n", grid_file);
    if(read_mesh(grid_file) == ERROR) {
      printf("octvce() - unfortunately read_mesh() generated an error.  Exiting ...\n");
      return(ERROR);
    }
  }
  else {
    /*start_time();*/
    refining_root_now = TRUE;
    if(refine(Root, num_times_refine_root, TRUE, TRUE, &Dummy_list) == NOERROR)
      refining_root_now = FALSE;
    else {
      printf("octvce() - Refining root cell unsuccessful.  Exiting ...\n");
      return(ERROR);
    }
    /*prn_time();*/
  }

  if((count_energy == TRUE) && (num_ICs > 0)) { /*Want to count total energies*/
    ec = malloc(sizeof(double)*num_ICs);
    vc = malloc(sizeof(double)*num_ICs);
    nc = malloc(sizeof(int)*num_ICs);
    
    for(i = 0; i < num_ICs; i++) {
      ec[i] = 0;
      nc[i] = 0;
      vc[i] = 0;
    }

    L = Leaves;
    while(L != NULL) {
      C = L -> cell_loc;

      if((C -> IC_flag != OUTSIDE) && (C -> IC_flag != UNKNOWN)) { /*So cell is inside some IC region*/
	i = (int) floor(C -> IC_flag);
	nc[i]++;
	if(twoD == 2) {
	  ec[i] += (C->Flow_data.State_vec.RhoE)*(C->centroid[1])*(C->cell_volume); /*Must get total J/rad*/
	  vc[i] += (C->centroid[1])*(C->cell_volume); /*Get total m^3/rad*/
	}
	else {
	  ec[i] += (C->Flow_data.State_vec.RhoE)*(C->cell_volume); /*Get total J*/
	  vc[i] += (C->cell_volume); /*Get total m^3*/
	}	
      }

      L = L ->next;
    }

    printf("\noctvce(): Here are the results from energy counting ...\n");
    
    for(i = 0; i < num_ICs; i++) {
      if(twoD == 2)
	printf("For IC region %d, tot no. cells is %d, tot energy is %9.8e J/rad, tot vol. is %9.8e m^3/rad\n", i, nc[i], ec[i], vc[i]);
      else if(twoD == 1)
	printf("For IC region %d, tot no. cells is %d, tot energy is %9.8e J/m, tot vol is %9.8e m^2\n", i, nc[i], ec[i], vc[i]);
      else printf("For IC region %d, tot no. cells is %d, tot energy is %9.8e J, tot vol is %9.8e m^3\n", i, nc[i], ec[i], vc[i]);
    }

    return(NOERROR);
  }

#if !NOTINCLUDE
  if(march_soln() != ERROR)
    printf("octvce() - Successfully marched the solution to its finality\n");      
  else {
    printf("octvce() - Couldn't march solution to its finality.  Something went wrong.  Exiting ...\n");
    return(ERROR);
  }
#endif

  close_history_files();

  return(NOERROR);
}

