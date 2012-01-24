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

/**\file Source file for solving Euler equations and computing related flow quantities for VCE octree cells*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ov_kernel.h"
#include "ov_flow.h"
#include "ov_adapt.h"
#include "ov_connect.h"
#include "ov_icbc.h"
#include "ov_io.h"
#include "ov_merge.h"
#include "ov_setgeom.h"
#include "ov_recon.h"
#include "ov_thermo.h"
#include "time_keeper.h"

#define NOTINCLUDE 0
#define DEBUG 0
#define GET_RID_EPS 0 /*Enforce want flow quantities whose absolute value < SMALL_NUM_LIMIT to go to 0*/
#define SMALL_NUM_LIMIT (DBL_EPSILON)
#define SOD_ON 0 /*If I'm running shock tube, I'll want to output solution at a given time and not at a given step*/
#define SOD_TIME 6.0e-4

extern short int twoD;
extern double twoD_thickness;
extern History_loc Histories;
extern Gas_dat gas_data;
extern Jwlb_coeff JWLB_coeffs;
extern State_vector IC_ambient;
extern Cfl_dat CFL_data;
extern short int flux_type;
extern double global_dt;
extern short int set_global_dt;
extern int number_of_threads;
extern int num_timesteps_run;
extern double finish_time;
extern double current_flow_time;
extern char many_limit;
extern char always_many_limit;
extern short int min_refinement_level;
extern short int max_refinement_level;
extern short int wall_rep;
extern double root_scale;
extern List_leaf Leaves;
extern List_leaf Leaves_tail;
extern List_leaf * Leaf_heads;
extern List_leaf * Leaf_tails;
extern double * min_dt;
extern int * num_cells_per_thread;
extern List_vtx_glob Vtxlist;
extern short int only_get_boundary_cells_flow;
extern short int want_to_adapt;
extern short int adapt_type;
extern int num_timesteps_adapt;
extern int num_timesteps_dump_data;
extern double time_interval_write;
extern double adapt_param;
extern double refine_param;
extern double coarsen_param;
extern double adapt_noise;
extern IC_heat ICheat;
extern Deton_pts Deton;
extern BC_data BC_est;
extern BC_data BC_wst;
extern BC_data BC_nth;
extern BC_data BC_sth;
extern BC_data BC_upr;
extern BC_data BC_lwr;

extern int stepnow;
extern int RK;

#if !NOTINCLUDE
extern int atime;
extern Flux_vector sflx;
List_leaf this_cell;
short int checking;
extern Cart_cell Root;
double flowtime;
#endif

#if SUPERSONIC_VORTEX
extern double M_in;
extern double max_drr_tol;
#endif

/*------------------------------------------------------------------*/

/**\brief March the flow soln - march through each timestep, compute and sum fluxes, output data and adapt
   if necessary.  Will march until num_timesteps_run reached*/

short int march_soln(void)
{
  int num_cells_flow;
  double current_time = 0;
  int run_timesteps_so_far = 0;
  int adapt_timesteps_so_far = 1;
  int data_dump_timesteps_so_far = 1;
  double write_time_so_far = 0;
  double next_time_write = time_interval_write;
  int history_dump_timesteps_so_far = 1;
  int num_refine = 0;
  int num_coarsen = 0;

#if NOTINCLUDE
  short int i, j;
  List_leaf L;
  Cart_cell C;
  int tmp;
#endif

#if !NOTINCLUDE
  stepnow = 0;
#endif

  /*Preliminaries at time 0*/
  if(run_timesteps_so_far == 0)
    { 
      /*Count the no. of plotting and computational cells in domain*/
      num_cells_flow = count_cells();
      
      /*Write the mesh and/or flow soln for time 0*/
#if VTXM
      if(write_soln(run_timesteps_so_far, num_cells_flow, 0) == ERROR)
	return(ERROR);
#endif
      
      if(output_history_loc(0) == ERROR)
	return(ERROR);

      /*Allocate flux vectors (do in serial as only done once for coarse grid)*/
      if(alloc_fluxes(Leaves, Leaves_tail, ALL) == ERROR)
	{
	  printf("march_soln(): Error allocating flux vectors\n");
	  return(ERROR);
	}

      /*Merge cells that need to be merged*/
      if(merge_all_cells() == ERROR)
	return(ERROR);

      /*Compute spatial gradients for next timestep*/
      get_grads_all_cells(Leaves, Leaves_tail, FALSE);

      /*Compute limiter values for each cell*/
      compute_limiter(Leaves, Leaves_tail, FALSE); 
    }

#if NOTINCLUDE
  tmp = num_timesteps_dump_data;
#endif

#if WARHEADBURN
  flowtime = 0;
#endif

  /*Now march in time*/
  while(((num_timesteps_run > 0) && (run_timesteps_so_far < num_timesteps_run)) || 
	((finish_time > 0) && (current_time < finish_time)))
    {       
      /*Compute global time step first*/
      if(set_global_dt == FALSE)
	get_global_timestep();
      
      /*Get flow solution for next timestep*/
#if NOTINCLUDE
      printf("Now printin' time 4 flow calc ...\n");
      start_time();
#endif

      if(get_flow_soln(&current_time) == ERROR)
	return(ERROR);
      /*CFL no. (hence global timestep) may be modified at this point*/

#if WARHEADBURN
     flowtime = current_time;
#endif

#if NOTINCLUDE
      prn_time();
#endif 

      /*Printout some useful info*/
      run_timesteps_so_far++;
      printf("Step %d, time %e s, dt %e, CFL %e, num_cells %d ", run_timesteps_so_far, current_time, global_dt, CFL_data.CFL,
	     num_cells_flow);      
      if((want_to_adapt == FALSE) || (adapt_timesteps_so_far != num_timesteps_adapt))
	printf("\n");
      
#if !NOTINCLUDE
      stepnow = run_timesteps_so_far;
#endif

#if SUPERSONIC_VORTEX
      if(calc_vortex_err(1.0, M_in, 1.0, run_timesteps_so_far, num_cells_flow, current_time) == TRUE)
	return(NOERROR);
#endif

      /*Commence adaptation if required - for the moment we can only refine ONCE each time i.e. 
	allow_multiple_refinements == FALSE, for simplicity's sake*/
#if NOTINCLUDE
      printf("Now printing time for ref ... \n");      
      start_time();
#endif

      if((want_to_adapt == TRUE) && (adapt_timesteps_so_far == num_timesteps_adapt)) 
	{  
	  if(adapt(&num_refine, &num_coarsen) == ERROR)
	    return(ERROR);

	  adapt_timesteps_so_far = 1; /*Reset*/
	  printf(", refined %d cells, coarsened %d cells\n", num_refine, num_coarsen);

	  if(Deton.deton_active == FALSE)
	    {
	      /*Merge cells that need to be merged after adapting*/
	      if(merge_all_cells() == ERROR)
		return(ERROR);
	  
	      /*Count the no. of plotting and computational cells in domain*/
	      num_cells_flow = count_cells(); /*Same no. cells as before (unless we just adapted)*/
	    }
	}
      else adapt_timesteps_so_far++;

      if(Deton.deton_active == TRUE) 
	{
	  /*Merge cells that need to be merged after adapting*/
	  if(merge_all_cells() == ERROR)
	    return(ERROR);
	  
	  /*Count the no. of plotting and computational cells in domain*/
	  num_cells_flow = count_cells(); /*Same no. cells as before (unless we just adapted)*/
	}

      /*Write solution data*/
      if(num_timesteps_dump_data > 0) {
	if(data_dump_timesteps_so_far == num_timesteps_dump_data) 
	  { 
	    if(write_soln(run_timesteps_so_far, num_cells_flow, current_time) == ERROR)
	      return(ERROR);
	    data_dump_timesteps_so_far = 1; /*Reset*/  
	  }
	else data_dump_timesteps_so_far++;
      }
      else { /*So we output by flow time, not flow timesteps*/
	write_time_so_far += global_dt; /*Always increase this*/

	if(write_time_so_far >= next_time_write) { 
	  if(write_soln(run_timesteps_so_far, num_cells_flow, current_time) == ERROR)
	    return(ERROR);
	  next_time_write += time_interval_write; 
	}
      }

      /*Write history files*/
      if(Histories.num_locs > 0)
	{
	  if(history_dump_timesteps_so_far == Histories.num_timesteps_history)
	    {
	      if(output_history_loc(current_time) == ERROR)
		return(ERROR);
	      history_dump_timesteps_so_far = 1; /*Reset*/
	    }
	  else history_dump_timesteps_so_far++;
	}

#if SOD_ON
      if(current_time > SOD_TIME) {
	num_cells_flow = count_cells();
	write_soln(run_timesteps_so_far, num_cells_flow, current_time);
	return(NOERROR);
      }
#endif

#if 0 
      if((stepnow == 100) || (stepnow == 200) || (stepnow == 300) || (stepnow == 400) || (stepnow == 500) || 
	 (stepnow == 600) || (stepnow == 700) || (stepnow == 800) || (stepnow == 900)) {
	num_cells_flow = count_cells();
	write_soln(run_timesteps_so_far, num_cells_flow, current_time);
      }      
#endif
            
#if NOTINCLUDE
      prn_time();
#endif 
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Helper function to count and/or flag the no. cells for different threads (in serial as very quick)*/

int count_cells(void)
{
  int num_cells = 0;
  List_leaf L;

  L = Leaves;
  while(L != NULL)
    {
#if DETON
      if(L -> cell_loc -> un_det != TRUE)
#endif 
	num_cells++;
      
      L = L -> next;
    }
  
  if(number_of_threads > 1)
    { 
      /*Get no. of threads for each cell (except last thread, which may have more/less than others)*/

      int i;
      int cells_per_thread = ROUND(((double) num_cells)/((double) number_of_threads)); 
      int cells = 1; /*Counter for no. cells in a sublist*/
      int thread_num = 0; 
      
      if(num_cells <= number_of_threads)
	{
	  /*Unlikely?  But reduce no. threads so num_cells == number_of_threads, then must allocate new sublist pointers*/
	}

      /*Now traverse list and tag cells as belonging to different sublists*/

      L = Leaves; 
      Leaf_heads[0] = L; 
      while(L != NULL)
	{
	  if(thread_num != (number_of_threads-1)) {
	    if(cells < cells_per_thread) {
	      L->thread_num = thread_num;

#if DETON
	      if(L->cell_loc->un_det != TRUE)
		cells++;
#endif
	    }
	    else if(cells == cells_per_thread) {
#if DETON
	      while((L != NULL) && (L->cell_loc->un_det == TRUE)) {
		L->thread_num=thread_num;
		L=L->next;
	      }
#endif
	      Leaf_tails[thread_num]=L;
	      Leaf_tails[thread_num]->thread_num = thread_num;

	      thread_num++; 
	      Leaf_heads[thread_num] = L->next; 
	      Leaf_heads[thread_num]->thread_num = thread_num;
	      L = Leaf_heads[thread_num];

#if DETON
	
	      if(L->cell_loc->un_det == TRUE) {
		while((L!=NULL) && (L->cell_loc->un_det == TRUE)) {
		  L->thread_num=thread_num;
		  L = L->next;
		}
		L = L->prev; /*Because at the very end of the while loop, L=L->next again*/
	      }
#endif
	      cells=2;
	    }
	  }
	  else if(thread_num == (number_of_threads-1)) {
	    L -> thread_num = thread_num;
	    if(L -> next == NULL)
	      Leaf_tails[thread_num] = L;
	  }
	    
	  L = L -> next;
	}

      for(i = 0; i < number_of_threads; i++)
	{
	  if(i == (number_of_threads - 1))
	    num_cells_per_thread[i] = num_cells - (number_of_threads-1)*cells_per_thread;
	  else num_cells_per_thread[i] = cells_per_thread;
	}
    }

  return(num_cells);
}

/*------------------------------------------------------------------*/

/**\brief Traverse list of leaf cells to get minimum timestep*/

void get_global_timestep(void)
{
  short int i;
  double mint;
  List_leaf node, Head, Tail;
  Cart_cell C;
  int thread_num = 0;

  for(i = 0; i < number_of_threads; i++)
    min_dt[i] = DBL_MAX; /*Set it to be impossibly large to begin with*/

  Head = Leaves;
  Tail = Leaves_tail;
  
  if(number_of_threads > 1)
    {
      #pragma omp parallel private(mint, node, Head, Tail, C, thread_num) 
      {	
        #ifdef _OPENMP
	thread_num = omp_get_thread_num();
        #endif

	Head = Leaf_heads[thread_num];
	Tail = Leaf_tails[thread_num];		

	node = Head;
	mint = min_dt[thread_num]; 

	while(1)
	  { 
	    C = node -> cell_loc;

	    if(
#if DETON
	       (C -> un_det != TRUE) &&
#endif
	       (C -> is_small_cell == FALSE))
	      {
		mint = min_local_cell_timestep(C);

		if(mint < min_dt[thread_num])
		  min_dt[thread_num] = mint;
	      }

	    if(node == Tail)
	      break;
	    else node = node -> next;	
	  }
      }
    }
  else
    {
      node = Head;
      mint = min_dt[0];

      while(1)
	{ 
	  C = node -> cell_loc;

	  if(
#if DETON
	     (C -> un_det != TRUE) &&
#endif
	     (C -> is_small_cell == FALSE))
	    {
	      mint = min_local_cell_timestep(C);
	      
	      if(mint < min_dt[0])
		min_dt[0] = mint;
	    }

	  if(node == Tail)
	    break;
	  else node = node -> next;	
	}
    }

  global_dt = min_dt[0]; /*Finally get the minimum out of all local minimums*/
  for(i = 0; i < number_of_threads; i++) 
    {
      if(min_dt[i] < global_dt)
	global_dt = min_dt[i];
    }
}

/*------------------------------------------------------------------*/

/**\brief Get the flow solution for the next time step (in serial/parallel)*/

short int get_flow_soln(double *time)
{  
  List_leaf Head, Tail; /*Head and tail of sublists processed by each thread*/
  int thread_num = 0;
  short int error_found = FALSE; /*Can't return in parallel construct*/  
  int i;
  BC_data *BC_now;
  BC_now = 0;
  if(ICheat.num_pts > 0)
    calc_trans_IC(*time, &ICheat, 1);
    
  for(i = 0; i < 6; i++) /*Get new transient border boundary values if necessary*/
    {
      switch(i)
	{
	case EST:
	  BC_now = &BC_est; break;
	case WST:
	  BC_now = &BC_wst; break;
	case NTH:
	  BC_now = &BC_nth; break;
	case STH:
	  BC_now = &BC_sth; break;
	case UPR:
	  BC_now = &BC_upr; break;
	case LWR:
	  BC_now = &BC_lwr; break;
	}

      if(BC_now -> type == TRANSIENT)
	calc_trans_BC(*time, BC_now, 1);	
    }

  Head = Leaves;
  Tail = Leaves_tail;
        
  if(number_of_threads > 1)
    {
      #pragma omp parallel private(Head, Tail, thread_num) 
      { 		
        #ifdef _OPENMP
	thread_num = omp_get_thread_num();
        #endif

	/*Now get sublists to work on*/
	Head = Leaf_heads[thread_num];
	Tail = Leaf_tails[thread_num];	
	RK = 1;

	/*Exchange fluxes and compute global time step*/
	if(all_cells_exchange_fluxes(Head, Tail, *time, 1) == ERROR) 
	  {
	    printf("get_flow_soln(): Error computing fluxes RK step 1\n");
	    error_found = TRUE;
            #pragma omp flush(error_found) /*Let error_found now be seen across all threads*/
	  }
	  
        #pragma omp barrier  /*Must ensure all threads compute fluxes before integrating in time to prevent race conditions*/

	/*Perform 1st stage of 2nd order Runge-Kutta method*/

	if(set_global_dt == TRUE)
	  {
	    if(error_found == FALSE)
	      {
		if(time_integrate(Head, Tail, 1, thread_num, FALSE) == ERROR)
		  {
		    printf("get_flow_soln(): Error time integrating RK step 1\n");
		    error_found = TRUE;
                    #pragma omp flush(error_found) 
		  }
	      }
	  }
	else 
	  {
	    if(error_found == FALSE)
	      {	
		checking = TRUE;
		CFL_data.CFL = CFL_data.CFL_max; /*First use the maximum CFL*/
		if(time_integrate(Head, Tail, 1, thread_num, TRUE) == ERROR)
		  {
		    printf("get_flow_soln(): Error time integrating RK step 1\n");
		    error_found = TRUE;
                    #pragma omp flush(error_found) 
		  }
		checking = FALSE; 
	      }

             #pragma omp barrier
	
	    if(error_found == FALSE)
	      {	
		if(CFL_data.CFL_ok == FALSE) /*CFL no. was revised*/
		  {
		    #pragma omp barrier

                    #pragma omp single
		    {
		      for(i = 0; i < number_of_threads; i++)
			if(CFL_data.CFL > CFL_data.CFL_min[i])
			  CFL_data.CFL = CFL_data.CFL_min[i];

#if 0
	printf("CFL_min - %e %e\n",CFL_data.CFL_min[0],CFL_data.CFL_min[1]);
#endif

		      CFL_data.CFL_ok = TRUE; /*By now the CFL no. should be acceptable*/
		    } 

                    #pragma omp barrier

		    if(time_integrate(Head, Tail, 1, thread_num, FALSE) == ERROR)
		      {
			printf("get_flow_soln(): Error time integrating RK step 1\n");
			error_found = TRUE;
                        #pragma omp flush(error_found) 
		      }
		  }
	      }	
	  }
	/*Do we need to compute the gradients and limiters again after each RK step????*/

        #pragma omp barrier

	/*Finally compute spatial gradients for next time step*/
	get_grads_all_cells(Head, Tail, FALSE);

	/*Finally decide if using multiple limiters or single limiter OK (if we only wanna use 1 limiter)*/
	if(always_many_limit == 'n') {
	  if(EQ(CFL_data.CFL, CFL_data.CFL_max))
	    many_limit = 'y';
	  else many_limit = 'n';
	}

	/*Compute limiter values for each cell*/
	compute_limiter(Head, Tail, FALSE);

        #pragma omp single
	{
	  if(ICheat.num_pts > 0)
	    calc_trans_IC((*time)+((CFL_data.CFL)*global_dt/2.0), &ICheat, 2);
	  
	  for(i = 0; i < 6; i++) /*Get new transient border boundary values if necessary*/
	    {
	      switch(i)
		{
		case EST:
		  BC_now = &BC_est; break;
		case WST:
		  BC_now = &BC_wst; break;
		case NTH:
		  BC_now = &BC_nth; break;
		case STH:
		  BC_now = &BC_sth; break;
		case UPR:
		  BC_now = &BC_upr; break;
		case LWR:
		  BC_now = &BC_lwr; break;
		}

	      if(BC_now -> type == TRANSIENT)
		calc_trans_BC((*time)+((CFL_data.CFL)*global_dt/2.0), BC_now, 2); /*For RK step 2*/		
	    }	  
	}
	  
        #pragma omp barrier 
	RK = 2;
	/*Exchange fluxes again with estimated new state*/
	if(error_found == FALSE)
	  {	
	    if(all_cells_exchange_fluxes(Head, Tail, *time, 2) == ERROR)
	      {
		printf("get_flow_soln(): Error computing fluxes RK step 2\n");
		error_found = TRUE;
                #pragma omp flush(error_found) 
	      }
	  }

        #pragma omp barrier 

	/*Perform 2nd stage of 2nd order Runge-Kutta method - time derivatives also computed here*/
	if(error_found == FALSE)
	  {	
	    if(time_integrate(Head, Tail, 2, thread_num, FALSE) == ERROR)
	      {
		printf("get_flow_soln(): Error time integrating RK step 2\n");
		error_found = TRUE;
                #pragma omp flush(error_found) 
	      }
	  }

        #pragma omp barrier

	/*Finally compute spatial gradients for next time step*/
	get_grads_all_cells(Head, Tail, FALSE);

	/*Compute limiter values for each cell*/
	compute_limiter(Head, Tail, FALSE);
      }
    }
  else
    {
      RK = 1;
      if(all_cells_exchange_fluxes(Head, Tail, *time, 1) == ERROR) 
	{
	  printf("get_flow_soln(): Error computing fluxes RK step 1\n");
	  return(ERROR);
	}

      if(set_global_dt == TRUE)
	{
	  if(time_integrate(Head, Tail, 1, 0, FALSE) == ERROR)
	    {
	      printf("get_flow_soln(): Error time integrating RK step 1\n");
	      return(ERROR);
	    } 
	}
      else 
	{
	  checking = TRUE;
	  CFL_data.CFL = CFL_data.CFL_max;
	  if(time_integrate(Head, Tail, 1, 0, TRUE) == ERROR)
	    {
	      printf("get_flow_soln(): Error time integrating RK step 1\n");
	      return(ERROR);
	    }
	  checking = FALSE;
	  
#if 1
	  if(CFL_data.CFL_ok == FALSE) /*CFL no. was revised*/
	    {	
	      CFL_data.CFL = CFL_data.CFL_min[0];
	      CFL_data.CFL_ok = TRUE; /*By now the CFL no. should be acceptable*/
	  
	      if(time_integrate(Head, Tail, 1, 0, FALSE) == ERROR)
		{
		  printf("get_flow_soln(): Error time integrating RK step 1\n");
		  return(ERROR);
		}
	    }
#else
	  CFL_data.CFL_ok = TRUE;
#endif
	}

      get_grads_all_cells(Head, Tail, FALSE);

      /*Finally decide if using multiple limiters or single limiter OK (if we only wanna use 1 limiter)*/
      if(always_many_limit == 'n') {
	if(EQ(CFL_data.CFL, CFL_data.CFL_max))
	  many_limit = 'y';
	else many_limit = 'n';
      }
      
      compute_limiter(Head, Tail, FALSE);            

      if(ICheat.num_pts > 0)
	calc_trans_IC((*time)+((CFL_data.CFL)*global_dt/2.0), &ICheat, 2);
	  
      for(i = 0; i < 6; i++) /*Get new transient border boundary values if necessary*/
	{
	  switch(i)
	    {
	    case EST:
	      BC_now = &BC_est; break;
	    case WST:
	      BC_now = &BC_wst; break;
	    case NTH:
	      BC_now = &BC_nth; break;
	    case STH:
	      BC_now = &BC_sth; break;
	    case UPR:
	      BC_now = &BC_upr; break;
	    case LWR:
	      BC_now = &BC_lwr; break;
	    }

	  if(BC_now -> type == TRANSIENT)
	    calc_trans_BC((*time)+((CFL_data.CFL)*global_dt/2.0), BC_now, 2); /*For RK step 2*/		
	}

      RK = 2;
      if(all_cells_exchange_fluxes(Head, Tail, *time, 2) == ERROR)
	{
	  printf("get_flow_soln(): Error computing fluxes RK step 2\n");
	  return(ERROR);
	}

      if(time_integrate(Head, Tail, 2, 0, FALSE) == ERROR)
	{
	  printf("get_flow_soln(): Error time integrating RK step 2\n");
	  return(ERROR);
	}
      
      get_grads_all_cells(Head, Tail, FALSE);
      compute_limiter(Head, Tail, FALSE);            
    }

  if(set_global_dt == FALSE)
    global_dt = global_dt*(CFL_data.CFL); /*Must alter global timestep so it incorporates CFL scaling factor*/
  *time = (*time) + global_dt; /*Advance to next timestep*/

  current_flow_time = *time;

#if DETON
  /*So this is where we check if the detonation wave(s) has now passed the undetonated cells*/
  if((Deton.num_deton_pts > 0) && (Deton.deton_active == TRUE)) {
    Deton.deton_active = FALSE;
    comp_new_deton_cells(*time, Leaves, Leaves_tail); /*If all cells have been detonated, deton_active remains FALSE*/
  }
#endif
  
#if CHECK_CRAP_QUANT
  Head = Leaves;
  while(Head != NULL) {

    if(Head->cell_loc->Flow_data.State_vec.Pres <= 0) {
      printf("-ve pressure encountered - %g\n", Head->cell_loc->Flow_data.State_vec.Pres);
      return(ERROR);
    }
    else if(Head->cell_loc->Flow_data.State_vec.Pres > 0) {
    }
    else {
      printf("Pressure %g encountered\n",Head->cell_loc->Flow_data.State_vec.Pres);
      return(ERROR);
    }

    if((Head->cell_loc->Flow_data.State_vec.RhoE - 0.5*(SQR(Head->cell_loc->Flow_data.State_vec.RhoU)+SQR(Head->cell_loc->Flow_data.State_vec.RhoV)+SQR(Head->cell_loc->Flow_data.State_vec.RhoW))/(Head->cell_loc->Flow_data.State_vec.Rhof)) <= 0) {
      printf("Energy low for cell %hd (%e %e %e), it's %e\n",Head->cell_loc->cell_level,Head->cell_loc->centroid[0],Head->cell_loc->centroid[1],Head->cell_loc->centroid[2],(Head->cell_loc->Flow_data.State_vec.RhoE - 0.5*(SQR(Head->cell_loc->Flow_data.State_vec.RhoU)+SQR(Head->cell_loc->Flow_data.State_vec.RhoV)+SQR(Head->cell_loc->Flow_data.State_vec.RhoW))/(Head->cell_loc->Flow_data.State_vec.Rhof)));
      printf("r %g,ru %g, rv %g, rw %g, re %g, rp %g\n",Head->cell_loc->Flow_data.State_vec.Rhof,Head->cell_loc->Flow_data.State_vec.RhoU,
	     Head->cell_loc->Flow_data.State_vec.RhoV,Head->cell_loc->Flow_data.State_vec.RhoW,Head->cell_loc->Flow_data.State_vec.RhoE,
	     Head->cell_loc->Flow_data.State_vec.rhof_prod);
      return(ERROR);
    }

    Head = Head -> next;
  }
#endif

  if(error_found == TRUE)
    return(ERROR);
  else return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Use Runge-Kutta method (for the moment only 2nd order) to sum fluxes then time integrate to get states at next timestep*/

short int time_integrate(List_leaf L, List_leaf Tail, short int stage, int thread_num, short int check_CFL)
{
  short int i, j, k, opp_dir, flux_sign;
  List_leaf node = L;
  List_merge merge_node;
  Cart_cell C, S;
  State_vector boundary_state;
  
  double sum_mss, sum_mss_prod;
  double sum_xmnt;
  double sum_ymnt;
  double sum_zmnt;
  double sum_e;
  double tot_area, vol, area, temp_pres, wall_norm[3];

#if AFTERBURN
  double sum_mss_prod_ub;
#endif

#if DETON
  double flux_area[4], flux_area_opp[4];
#endif

  double max_dr = 0;
  double max_dp = 0; /*Maximum relative changes in density and pressure (for revision of CFL no.)*/
  double dr, dp;

#if NOTINCLUDE
  int sc_count;
  Cart_cell temp_neighb;
  int num_smalls;
  short int found = 0;
  double tmp1, tmp2;
#endif

  double r_val, r_L, r_R;
  short int is_corner_cell;
  area = 0; temp_pres = 0; r_val = 0;
  while(node != NULL)
    {
      C = node -> cell_loc;

      if(
#if DETON
	 (C -> un_det != TRUE) &&
#endif
	 (C -> is_small_cell == FALSE)) /*Only deal with 'normal' cells; normal cells part of a merged cluster will
					  sum fluxes over all linked cells*/
	{
	  if((check_CFL == FALSE) && (stage == 1) && (set_global_dt == FALSE)) /*This happens when CFL no. is revised, else check_CFL is 
										 TRUE.  The flux sums have already been computed before*/
            {
	      sum_mss = C -> Flow_data.Flx_sums.sum_mss;
	      sum_xmnt = C -> Flow_data.Flx_sums.sum_xmnt;
	      sum_ymnt = C -> Flow_data.Flx_sums.sum_ymnt;
	      sum_zmnt = C -> Flow_data.Flx_sums.sum_zmnt;
	      sum_e = C -> Flow_data.Flx_sums.sum_e;
	      sum_mss_prod = C -> Flow_data.Flx_sums.sum_mss_prod;
#if AFTERBURN
	      sum_mss_prod_ub = C -> Flow_data.Flx_sums.sum_mss_prod_ub;
#endif
	    }
	  else
	    {
	      sum_mss = 0; 
	      sum_xmnt = 0;
	      sum_ymnt = 0;
	      sum_zmnt = 0;
	      sum_e = 0;
	      sum_mss_prod = 0;
#if AFTERBURN
	      sum_mss_prod_ub = 0;
#endif

	      /*Sum cell to cell fluxes*/
	      for(i = 0; i < NUM_FACES; i++) 
		{
		  if((i == EST) || (i == NTH) || (i == UPR) || (C -> face_neighbours[i][0] == NULL)) /*Sign of flux vectors*/
		    flux_sign = 1;
		  else flux_sign = -1;	      

		  if(((C -> face_neighbours[i][0] == NULL) || 
		      ((C -> face_neighbours[i][1] == NULL) 
#if DETON
		       && (C -> face_neighbours[i][0] -> un_det != TRUE) 
#endif
		       )) &&
		     (C -> Flow_data.Face_fluxes[i][0] != NULL))
		    { /*1 valid flux vector in this direction (won't exist between linked cells) - use whole face*/
		  
		      tot_area = (C -> flux_area[i][0]) + (C -> flux_area[i][1]) + (C -> flux_area[i][2]) + (C -> flux_area[i][3]);

		      if(twoD == 2) {
			if(i == NTH)
			  r_val = (C->verticies[1]->loc[1]);
			else if(i == STH)
			  r_val = C->verticies[0]->loc[1];
			else {
#if AXI_OK
			  r_val = ((C->flux_area[i][0])*(C -> face_r_centroid[i][0])+(C->flux_area[i][1])*(C -> face_r_centroid[i][1]))/
			    (C->flux_area[i][0] + C->flux_area[i][1]);		
#else
			  r_val = C->centroid[1];
#endif
			}
		      }
		      else r_val = 1;

		      sum_mss += (C -> Flow_data.Face_fluxes[i][0] -> Mss_flx)*tot_area*flux_sign*r_val;
		      sum_xmnt += (C -> Flow_data.Face_fluxes[i][0] -> X_mntm_flx)*tot_area*flux_sign*r_val;
		      sum_ymnt += (C -> Flow_data.Face_fluxes[i][0] -> Y_mntm_flx)*tot_area*flux_sign*r_val;
		      sum_zmnt += (C -> Flow_data.Face_fluxes[i][0] -> Z_mntm_flx)*tot_area*flux_sign*r_val;
		      sum_e += (C -> Flow_data.Face_fluxes[i][0] -> E_flx)*tot_area*flux_sign*r_val;

		      if(gas_data.num_species > 1) {
			sum_mss_prod += (C -> Flow_data.Face_fluxes[i][0] -> mss_flx_prod)*tot_area*flux_sign*r_val;
#if AFTERBURN
			sum_mss_prod_ub += (C->Flow_data.Face_fluxes[i][0] -> mss_flx_prod_ub)*tot_area*flux_sign*r_val;
#endif
		      }

#if 0
		      if(sum_mss <= 0) {
		      }
		      else if(sum_mss > 0) {
		      }
		      else {
			printf("1 neighb, Sum_mss is %g, dir %hd, C level %hd (%e %e %e), r_val is %g, flx is %g, area %g\n", sum_mss,i,C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],
			       r_val,C->Flow_data.Face_fluxes[i][0]->Mss_flx,tot_area);
			abort();
		      }
#endif
		    }
		  else if(C -> face_neighbours[i][1] != NULL) /*May have 4 flux vectors here*/
		    {
		      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
			{
			  if(
#if DETON
			     (C -> face_neighbours[i][j] -> un_det != TRUE) &&
#endif
			     (C -> Flow_data.Face_fluxes[i][j] != NULL))
			    { /*1 valid flux vector in this quadrant*/

			      if(twoD == 2) {
				if(i == NTH)
				  r_val = (C->verticies[1]->loc[1]);
				else if(i == STH)
				  r_val = C->verticies[0]->loc[1];
				else {
#if AXI_OK
				  if(j == 0)
				    r_val = C -> face_r_centroid[i][0];
				  else if(j == 1)
				    r_val = C -> face_r_centroid[i][1];		
#else
				  if(j == 0)
				    r_val = ((C->centroid[1])+(C->verticies[0]->loc[1]))/2.0;
				  else if(j == 1)
				    r_val = ((C->centroid[1])+(C->verticies[1]->loc[1]))/2.0;
#endif
				}
			      }
			      else r_val = 1;

			      sum_mss += (C -> Flow_data.Face_fluxes[i][j] -> Mss_flx)*(C -> flux_area[i][j])*flux_sign*r_val;
			      sum_xmnt += (C -> Flow_data.Face_fluxes[i][j] -> X_mntm_flx)*(C -> flux_area[i][j])*flux_sign*r_val;
			      sum_ymnt += (C -> Flow_data.Face_fluxes[i][j] -> Y_mntm_flx)*(C -> flux_area[i][j])*flux_sign*r_val;
			      sum_zmnt += (C -> Flow_data.Face_fluxes[i][j] -> Z_mntm_flx)*(C -> flux_area[i][j])*flux_sign*r_val;
			      sum_e += (C -> Flow_data.Face_fluxes[i][j] -> E_flx)*(C -> flux_area[i][j])*flux_sign*r_val;

			      if(gas_data.num_species > 1) {
				sum_mss_prod += (C -> Flow_data.Face_fluxes[i][j] -> mss_flx_prod)*(C -> flux_area[i][j])*flux_sign*r_val;
#if AFTERBURN
				sum_mss_prod_ub += (C->Flow_data.Face_fluxes[i][j]->mss_flx_prod_ub)*(C -> flux_area[i][j])*
				  flux_sign*r_val;
#endif
			      }
			    }
			}
		    }
		}
	  
	      /*Wall fluxes use special Obstructed_flux vector per obstructed cell*/
	      /*Even if 'smooth' surface is used, there might be fluid cells with wall norms that are NOT axis unit vectors i.e.
		mostly grid-aligned corner cells.  It would then be better to treat them with staircased representation.*/

	      is_corner_cell = FALSE;
	      if((twoD == 2) && (C -> cell_type == FLUID) && (((ABS(C->wall_norm[0]) < 1) && (ABS(C->wall_norm[0]) > AREA_DIFF_TOL)) && 
							      ((ABS(C->wall_norm[1]) < 1) && (ABS(C->wall_norm[1]) > AREA_DIFF_TOL)))) 
		is_corner_cell = TRUE;
	      	      
	      if((wall_rep == 1) && (is_corner_cell == FALSE))
		{
		  if(C -> Flow_data.Obstructed_flux != NULL)
		    {
		      if(twoD == 2) {
#if AXI_OK
			/*First get average radial coordinate of east and west unobstructed side lengths*/
			if((C->flux_area[WST][0] + C->flux_area[WST][1]) > 0) {
			  r_L = ((C->flux_area[WST][0])*(C->face_r_centroid[WST][0])+
				 (C->flux_area[WST][1])*(C->face_r_centroid[WST][1]))/(C->flux_area[WST][0] + C->flux_area[WST][1]);
			  r_L = r_L - (C -> verticies[0] -> loc[1]); /*Shift reference frame to lowest cell corner*/
			}
			else r_L = 0;

			if((C->flux_area[EST][0] + C->flux_area[EST][1]) > 0) {
			  r_R = ((C->flux_area[EST][0])*(C->face_r_centroid[EST][0])+
			     (C->flux_area[EST][1])*(C->face_r_centroid[EST][1]))/(C->flux_area[EST][0] + C->flux_area[EST][1]);
			  r_R = r_R - (C -> verticies[0] -> loc[1]);
			}
			else r_R = 0;

			if(C -> wall_norm[1] >= 0)  /*If wall normal is upward, then it's very simple*/
			  r_val = r_L + r_R + (C -> verticies[0] -> loc[1]);
			else { /*If wall normal is downward, need to apply transformation r' = l - r*/
			  if(r_L > 0)
			    r_L = (C -> cell_length) - r_L;
			  if(r_R > 0)
			    r_R = (C -> cell_length) - r_R;
			  /*If either r_L or r_R is 0, don't need to transform*/

			  r_val = (C -> cell_length) - (r_L + r_R) + (C -> verticies[0] -> loc[1]); /*Transform back*/
			}
#else
			if(
#if FLOAT_EQ_STABLE
			   EQ((C->flux_area[STH][0]+C->flux_area[STH][1]+C->flux_area[STH][2]+C->flux_area[STH][3]), 0)
#else
			   (C->flux_area[STH][0]+C->flux_area[STH][1]+C->flux_area[STH][2]+C->flux_area[STH][3]) == 0
#endif
			   )
			  r_val = C->verticies[0]->loc[1];
			else if(
#if FLOAT_EQ_STABLE
				EQ((C->flux_area[NTH][0]+C->flux_area[NTH][1]+C->flux_area[NTH][2]+C->flux_area[NTH][3]) , 0)
#else
				(C->flux_area[NTH][0]+C->flux_area[NTH][1]+C->flux_area[NTH][2]+C->flux_area[NTH][3]) == 0
#endif
				)
			  r_val = C->verticies[1]->loc[1];
			else r_val = C->centroid[1];
#endif
		      }
		      else r_val = 1;

		      sum_mss += (C -> Flow_data.Obstructed_flux -> Mss_flx)*(C->wall_area)*r_val;
		      sum_xmnt += (C -> Flow_data.Obstructed_flux -> X_mntm_flx)*(C->wall_area)*r_val;
		      sum_ymnt += (C -> Flow_data.Obstructed_flux -> Y_mntm_flx)*(C->wall_area)*r_val;
		      sum_zmnt += (C -> Flow_data.Obstructed_flux -> Z_mntm_flx)*(C->wall_area)*r_val;
		      sum_e += (C -> Flow_data.Obstructed_flux -> E_flx)*(C->wall_area)*r_val;

		      if(gas_data.num_species > 1) {
			sum_mss_prod += (C -> Flow_data.Obstructed_flux -> mss_flx_prod)*(C->wall_area)*r_val;
#if AFTERBURN
			sum_mss_prod_ub += (C -> Flow_data.Obstructed_flux -> mss_flx_prod_ub)*(C->wall_area)*r_val;
#endif
		      }
		    }
		}
	      else /*Choose a 'staircased' surface representation*/
		{		  
		  for(j = 0; j < 3; j++)
		    {
		      i = j*2; /*'Positive' directions along the 3 axes*/
		      opp_dir = i+1; /*'Negative' directions along the 3 axes*/

#if DETON
		      if(Deton.num_deton_pts > 0) {
			for(k = 0; k < 4; k++) {
			  flux_area[k] = C->flux_area[i][k];
			  flux_area_opp[k] = C->flux_area[opp_dir][k];
			}

			if(C -> face_neighbours[i][0] != NULL) {
			  if((C -> face_neighbours[i][1] == NULL)) {
			    if(C -> face_neighbours[i][0] -> un_det == TRUE) {		      
			      for(k = 0; k < 4; k++)
				flux_area[k] = 0; /*Undetonated cell is like a solid cell*/
			    }
			  }
			  else {
			    for(k = 0; k < 4; k++) {
			      if(C -> face_neighbours[i][k] -> un_det == TRUE)
				flux_area[k] = 0;
			    }
			  }
			}

			if(C -> face_neighbours[opp_dir][0] != NULL) {
			  if((C -> face_neighbours[opp_dir][1] == NULL)) {
			    if(C -> face_neighbours[opp_dir][0] -> un_det == TRUE) {		      
			      for(k = 0; k < 4; k++)
				flux_area_opp[k] = 0; 
			    }
			  }
			  else {
			    for(k = 0; k < 4; k++) {
			      if(C -> face_neighbours[opp_dir][k] -> un_det == TRUE)
				flux_area_opp[k] = 0;
			    }
			  }
			}
			
			tot_area = (flux_area_opp[0]+flux_area_opp[1]+flux_area_opp[2]+flux_area_opp[3]) -
			  (flux_area[0]+flux_area[1]+flux_area[2]+flux_area[3]);
		      }
		      else {
			tot_area = (C->flux_area[opp_dir][0]+C->flux_area[opp_dir][1]+C->flux_area[opp_dir][2]+C->flux_area[opp_dir][3]) -
			  (C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]);
		      }
#else
		      tot_area = (C->flux_area[opp_dir][0]+C->flux_area[opp_dir][1]+C->flux_area[opp_dir][2]+C->flux_area[opp_dir][3]) -
			(C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]);
		      /*Net obstructed area in the direction*/
#endif

		      if(ABS(tot_area) > 0) /*So a net obstruction does exist*/
			{
			  if(C -> Flow_data.Obstructed_flux == NULL) /*Allocate if not allocated yet*/
			    {
			      C -> Flow_data.Obstructed_flux = malloc(sizeof(struct flux_vector));
			      if(C -> Flow_data.Obstructed_flux == NULL)
				return(ERROR);
			    }

			  wall_norm[0] = 0; wall_norm[1] = 0; wall_norm[2] = 0;

			  if(tot_area < 0) /*Negative face has greater obstruction*/
			    wall_norm[j] = -1;
			  else wall_norm[j] = 1;

			  /*Fixer for 'corner cells'*/
			  r_val = 1;
			  if(twoD == 2) {
			    if(j == 0) 
			      r_val = C->centroid[1];
			    else if(j == 1) {
			      if(wall_norm[j] > 0)
				r_val = C->verticies[7]->loc[1]; /*Obstruction on north face*/
			      else r_val = C->verticies[0]->loc[1]; /*Obstruction on south face*/
			    }			    
			  }

			  boundary_state = set_boundary_state(&(C->Flow_data.State_vec), NULL, wall_norm, i, WALL, 
							      IRRELEVANT, IRRELEVANT, IRRELEVANT, NULL, 0);

			  if(flux_type == ADAPTIVE) {
			    if(C -> shocked_cell == 'y')
			      compute_fluxes(&(C->Flow_data.State_vec),&boundary_state,wall_norm,C->Flow_data.Obstructed_flux,EFM);
			    else compute_fluxes(&(C->Flow_data.State_vec),&boundary_state,wall_norm,C->Flow_data.Obstructed_flux,AUSMDV);
			  }
			  else compute_fluxes(&(C->Flow_data.State_vec),&boundary_state,wall_norm,C->Flow_data.Obstructed_flux,flux_type);

			  tot_area = ABS(tot_area);		  
			  sum_mss += (C -> Flow_data.Obstructed_flux -> Mss_flx)*tot_area*r_val;
			  sum_xmnt += (C -> Flow_data.Obstructed_flux -> X_mntm_flx)*tot_area*r_val;
			  sum_ymnt += (C -> Flow_data.Obstructed_flux -> Y_mntm_flx)*tot_area*r_val;
			  sum_zmnt += (C -> Flow_data.Obstructed_flux -> Z_mntm_flx)*tot_area*r_val;
			  sum_e += (C -> Flow_data.Obstructed_flux -> E_flx)*tot_area*r_val;

#if 0
			  if(EQ(C->centroid[0],0.0625) && EQ(C->centroid[1],0.0625)) {
			    printf("dir %hd, rv %g, ta %g, sm %g, sx %g, sy %g, sz %g, se %g\n",j,r_val,tot_area,
				   (C -> Flow_data.Obstructed_flux -> Mss_flx),
				   (C -> Flow_data.Obstructed_flux -> X_mntm_flx),
				   (C -> Flow_data.Obstructed_flux -> Y_mntm_flx),
				   (C -> Flow_data.Obstructed_flux -> Z_mntm_flx),
				   (C -> Flow_data.Obstructed_flux -> E_flx));
			  }
#endif

			  if(gas_data.num_species > 1) {
			    sum_mss_prod += (C -> Flow_data.Obstructed_flux -> mss_flx_prod)*tot_area*r_val;
#if AFTERBURN
			    sum_mss_prod_ub += (C -> Flow_data.Obstructed_flux -> mss_flx_prod_ub)*tot_area*r_val;
#endif
			  }
			}
		    }
		}
	     
	      if(C -> Merge_data != NULL) /*Then sum fluxes across all interfaces; normal cell has already been dealt with*/
		{
		  merge_node = C -> Merge_data -> Linked_cells;
	      
		  while(merge_node != NULL)
		    {
		      S = merge_node -> cell_loc;

		      if(S -> is_small_cell == TRUE) /*Sum fluxes*/
			{
			  /*Cell to cell fluxes*/
			  for(i = 0; i < NUM_FACES; i++) 
			    {
			      if((i == EST) || (i == NTH) || (i == UPR) || (S -> face_neighbours[i][0] == NULL))
				flux_sign = 1; /*Flux 'sign' always positive at border faces*/
			      else flux_sign = -1;

			      if(((S -> face_neighbours[i][0] == NULL) || 
				  ((S -> face_neighbours[i][1] == NULL)
#if DETON
				   && (S -> face_neighbours[i][0] -> un_det != TRUE)
#endif
				   )) && 
				 (S -> Flow_data.Face_fluxes[i][0] != NULL))
				{ 	
				  tot_area = (S -> flux_area[i][0]) + (S -> flux_area[i][1]) + (S -> flux_area[i][2]) + 
				    (S -> flux_area[i][3]);
		  
				  if(twoD == 2) {
				    if(i == NTH)
				      r_val = (S->verticies[1]->loc[1]);
				    else if(i == STH)
				      r_val = S->verticies[0]->loc[1];
				    else {
#if AXI_OK
				      r_val = ((S->flux_area[i][0])*(S -> face_r_centroid[i][0])+
					       (S->flux_area[i][1])*(S -> face_r_centroid[i][1]))/
					(S->flux_area[i][0] + S->flux_area[i][1]);	
#else
				      r_val = S->centroid[1];
#endif
				    }
				  }
				  else r_val = 1;

				  sum_mss += (S -> Flow_data.Face_fluxes[i][0] -> Mss_flx)*tot_area*flux_sign*r_val;
				  sum_xmnt += (S -> Flow_data.Face_fluxes[i][0] -> X_mntm_flx)*tot_area*flux_sign*r_val;
				  sum_ymnt += (S -> Flow_data.Face_fluxes[i][0] -> Y_mntm_flx)*tot_area*flux_sign*r_val;
				  sum_zmnt += (S -> Flow_data.Face_fluxes[i][0] -> Z_mntm_flx)*tot_area*flux_sign*r_val;
				  sum_e += (S -> Flow_data.Face_fluxes[i][0] -> E_flx)*tot_area*flux_sign*r_val;

				  if(gas_data.num_species > 1) {
				    sum_mss_prod += (S -> Flow_data.Face_fluxes[i][0] -> mss_flx_prod)*tot_area*flux_sign*r_val;
#if AFTERBURN
				    sum_mss_prod_ub += (S -> Flow_data.Face_fluxes[i][0] -> mss_flx_prod_ub)*tot_area*flux_sign*r_val;
#endif
				  }
				}
			      else if(S -> face_neighbours[i][1] != NULL) 
				{
				  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
				    {
				      if(
#if DETON
					 (S -> face_neighbours[i][j] -> un_det != TRUE) &&
#endif
					 (S -> Flow_data.Face_fluxes[i][j] != NULL))
					{ 
					  
					  if(twoD == 2) {
					    if(i == NTH)
					      r_val = (S->verticies[1]->loc[1]);
					    else if(i == STH)
					      r_val = S->verticies[0]->loc[1];
					    else {
#if AXI_OK
					      if(j == 0)
						r_val = S -> face_r_centroid[i][0];
					      else if(j == 1)
						r_val = S -> face_r_centroid[i][1];
#else
					      if(j == 0)
						r_val = ((S->centroid[1])+(S->verticies[0]->loc[1]))/2.0;
					      else if(j == 1)
						r_val = ((S->centroid[1])+(S->verticies[1]->loc[1]))/2.0;
#endif
					    }
					  }
					  else r_val = 1;

					  sum_mss += (S -> Flow_data.Face_fluxes[i][j] -> Mss_flx)*(S -> flux_area[i][j])*flux_sign*r_val;
					  sum_xmnt += (S -> Flow_data.Face_fluxes[i][j] -> X_mntm_flx)*(S -> flux_area[i][j])*
					    flux_sign*r_val;
					  sum_ymnt += (S -> Flow_data.Face_fluxes[i][j] -> Y_mntm_flx)*(S -> flux_area[i][j])*
					    flux_sign*r_val;
					  sum_zmnt += (S -> Flow_data.Face_fluxes[i][j] -> Z_mntm_flx)*(S -> flux_area[i][j])*
					    flux_sign*r_val;
					  sum_e += (S -> Flow_data.Face_fluxes[i][j] -> E_flx)*(S -> flux_area[i][j])*flux_sign*r_val;

					  if(gas_data.num_species > 1) {
					    sum_mss_prod += (S->Flow_data.Face_fluxes[i][j]->mss_flx_prod)*(S->flux_area[i][j])*
					      flux_sign*r_val;
#if AFTERBURN
					    sum_mss_prod_ub += (S->Flow_data.Face_fluxes[i][j]->mss_flx_prod_ub)*
					      (S->flux_area[i][j])*flux_sign*r_val;
#endif
					  }
					}
				    }
				}
			    }
		      
			  /*Cell to wall fluxes*/

			  is_corner_cell = FALSE;
			  if((twoD == 2) && (S -> cell_type == FLUID) && 
			     (((ABS(S->wall_norm[0]) < 1) && (ABS(S->wall_norm[0]) > AREA_DIFF_TOL)) && 
			      ((ABS(S->wall_norm[1]) < 1) && (ABS(S->wall_norm[1]) > AREA_DIFF_TOL))))
			    is_corner_cell = TRUE;
			  
			  if((wall_rep == 1) && (is_corner_cell == FALSE))
			    {
			      if(S -> Flow_data.Obstructed_flux != NULL)
				{
				  if(twoD == 2) {
#if AXI_OK
				    if((S->flux_area[WST][0] + S->flux_area[WST][1]) > 0) {
				      r_L = 
					((S->flux_area[WST][0])*(S->face_r_centroid[WST][0])+
					 (S->flux_area[WST][1])*(S->face_r_centroid[WST][1]))/
					(S->flux_area[WST][0] + S->flux_area[WST][1]);
				      r_L = r_L - (S -> verticies[0] -> loc[1]); /*Shift reference frame to lowest cell corner*/
				    }
				    else r_L = 0;

				    if((S->flux_area[EST][0] + S->flux_area[EST][1]) > 0) {
				      r_R = 
					((S->flux_area[EST][0])*(S->face_r_centroid[EST][0])+
					 (S->flux_area[EST][1])*(S->face_r_centroid[EST][1]))/
					(S->flux_area[EST][0] + S->flux_area[EST][1]);
				      r_R = r_R - (S -> verticies[0] -> loc[1]);
				    }
				    else r_R = 0;

				    if((S->flux_area[EST][0] + S->flux_area[EST][1]) > 0) {
				      r_R = ((S->flux_area[EST][0])*(S->face_r_centroid[EST][0])+
					     (S->flux_area[EST][1])*(S->face_r_centroid[EST][1]))/
					(S->flux_area[EST][0] + S->flux_area[EST][1]);
				      r_R = r_R - (S -> verticies[0] -> loc[1]);
				    }
				    else r_R = 0;

				    if(S -> wall_norm[1] >= 0)  
				      r_val = r_L + r_R + (S -> verticies[0] -> loc[1]);
				    else { 
				      if(r_L > 0)
					r_L = (S -> cell_length) - r_L;
				      if(r_R > 0)
					r_R = (S -> cell_length) - r_R;
				     				      
				      r_val = (S -> cell_length) - (r_L + r_R) + (S -> verticies[0] -> loc[1]); /*Transform back*/
				    }
#else
				    if((S->flux_area[STH][0]+S->flux_area[STH][1]+S->flux_area[STH][2]+S->flux_area[STH][3]) == 0)
				      r_val = S->verticies[0]->loc[1];
				    else if((S->flux_area[NTH][0]+S->flux_area[NTH][1]+S->flux_area[NTH][2]+S->flux_area[NTH][3]) == 0)
				      r_val = S->verticies[1]->loc[1];
				    else r_val = S->centroid[1];
#endif
				  }
				  else r_val = 1;

				  sum_mss += (S -> Flow_data.Obstructed_flux -> Mss_flx)*(S->wall_area)*r_val;
				  sum_xmnt += (S -> Flow_data.Obstructed_flux -> X_mntm_flx)*(S->wall_area)*r_val;
				  sum_ymnt += (S -> Flow_data.Obstructed_flux -> Y_mntm_flx)*(S->wall_area)*r_val;
				  sum_zmnt += (S -> Flow_data.Obstructed_flux -> Z_mntm_flx)*(S->wall_area)*r_val;
				  sum_e += (S -> Flow_data.Obstructed_flux -> E_flx)*(S->wall_area)*r_val;

				  if(gas_data.num_species > 1) {
				    sum_mss_prod += (S -> Flow_data.Obstructed_flux -> mss_flx_prod)*(S->wall_area)*r_val;
#if AFTERBURN
				    sum_mss_prod_ub += (S -> Flow_data.Obstructed_flux -> mss_flx_prod_ub)*(S->wall_area)*r_val;
#endif
				  }
				}
			    }
			  else 
			    {
			      for(j = 0; j < 3; j++)
				{
				  i = j*2;
				  opp_dir = i+1;

#if DETON
				  if(Deton.num_deton_pts > 0) {
				    for(k = 0; k < 4; k++) { /*First get the 'proper' flux areas*/
				      flux_area[k] = S->flux_area[i][k];
				      flux_area_opp[k] = S->flux_area[opp_dir][k];
				    }

				    /*Now change the flux areas accordingly*/
				    if(S -> face_neighbours[i][0] != NULL) {
				      if((S -> face_neighbours[i][1] == NULL)) {
					if(S -> face_neighbours[i][0] -> un_det == TRUE) {		      
					  for(k = 0; k < 4; k++)
					    flux_area[k] = 0; /*Undetonated cell is like a solid cell*/
					}
				      }
				      else {
					for(k = 0; k < 4; k++) {
					  if(S -> face_neighbours[i][k] -> un_det == TRUE)
					    flux_area[k] = 0;
					}
				      }
				    }

				    if(S -> face_neighbours[opp_dir][0] != NULL) {
				      if((S -> face_neighbours[opp_dir][1] == NULL)) {
					if(S -> face_neighbours[opp_dir][0] -> un_det == TRUE) {		      
					  for(k = 0; k < 4; k++)
					    flux_area_opp[k] = 0; 
					}
				      }
				      else {
					for(k = 0; k < 4; k++) {
					  if(S -> face_neighbours[opp_dir][k] -> un_det == TRUE)
					    flux_area_opp[k] = 0;
					}
				      }
				    }

				    tot_area = (flux_area_opp[0]+flux_area_opp[1]+flux_area_opp[2]+flux_area_opp[3]) -
				      (flux_area[0]+flux_area[1]+flux_area[2]+flux_area[3]);
				  }
				  else {
				    tot_area = (S->flux_area[opp_dir][0]+S->flux_area[opp_dir][1]+S->flux_area[opp_dir][2]+
						S->flux_area[opp_dir][3])
				      - (S->flux_area[i][0] + S->flux_area[i][1] + S->flux_area[i][2] + S->flux_area[i][3]);
				  }
#else
				  tot_area = (S->flux_area[opp_dir][0]+S->flux_area[opp_dir][1]+S->flux_area[opp_dir][2]+
					      S->flux_area[opp_dir][3])
				    - (S->flux_area[i][0] + S->flux_area[i][1] + S->flux_area[i][2] + S->flux_area[i][3]);
#endif

				  if(ABS(tot_area) > 0)
				    {
				      if(S -> Flow_data.Obstructed_flux == NULL)
					{
					  S -> Flow_data.Obstructed_flux = malloc(sizeof(struct flux_vector));
					  if(S -> Flow_data.Obstructed_flux == NULL)
					    return(ERROR);
					}

				      wall_norm[0] = 0; wall_norm[1] = 0; wall_norm[2] = 0;
				      if(tot_area < 0) /*Negative face has greater obstruction*/
					wall_norm[j] = -1;
				      else wall_norm[j] = 1;

				      /*Fixer for 'corner cells'*/
				      r_val = 1;
				      if(twoD == 2) {
					if(j == 0) 
					  r_val = S->centroid[1];
					else if(j == 1) {
					  if(wall_norm[j] > 0)
					    r_val = S->verticies[7]->loc[1]; /*Obstruction on north face*/
					  else r_val = S->verticies[0]->loc[1]; /*Obstruction on south face*/
					}			    
				      }

				      boundary_state = set_boundary_state(&(S->Flow_data.State_vec),NULL,wall_norm,i,WALL,
									  IRRELEVANT, IRRELEVANT, IRRELEVANT, NULL, 0);

				      if(flux_type == ADAPTIVE) {
					if(C -> shocked_cell == 'y')
					  compute_fluxes(&(S->Flow_data.State_vec), &boundary_state, wall_norm, 
							 S->Flow_data.Obstructed_flux, EFM);
					else compute_fluxes(&(S->Flow_data.State_vec), &boundary_state, wall_norm, 
							  S->Flow_data.Obstructed_flux, AUSMDV);
				      }
				      else compute_fluxes(&(S->Flow_data.State_vec), &boundary_state, wall_norm, 
							  S->Flow_data.Obstructed_flux, flux_type);

				      tot_area = ABS(tot_area);
				      sum_mss += (S -> Flow_data.Obstructed_flux -> Mss_flx)*tot_area*r_val;
				      sum_xmnt += (S -> Flow_data.Obstructed_flux -> X_mntm_flx)*tot_area*r_val;
				      sum_ymnt += (S -> Flow_data.Obstructed_flux -> Y_mntm_flx)*tot_area*r_val;
				      sum_zmnt += (S -> Flow_data.Obstructed_flux -> Z_mntm_flx)*tot_area*r_val;
				      sum_e += (S -> Flow_data.Obstructed_flux -> E_flx)*tot_area*r_val;

				      if(gas_data.num_species > 1) {
					sum_mss_prod += (S -> Flow_data.Obstructed_flux -> mss_flx_prod)*tot_area*r_val;
#if AFTERBURN
					sum_mss_prod_ub += (S -> Flow_data.Obstructed_flux -> mss_flx_prod_ub)*tot_area*r_val;
#endif
				      }
				    }
				}
			    }
			} /*Normal cell already visited*/
		  
		      merge_node = merge_node -> next;
		    }

#if NOTINCLUDE
		  if(num_smalls > 2)
		    printf("More than 2 cells in this cluster\n");
#endif
		}
      	    }

	  if(C -> Merge_data == NULL)
	    vol = C -> cell_volume;
	  else vol = C -> Merge_data -> cell_volume;

	  if(twoD != 0) { /*Extra precautionary measures*/
	    C -> Flow_data.State_vec.RhoW = 0;
	    sum_zmnt = 0; 
	  }

	  if(stage == 1) /*First stage of Runge-Kutta method*/
	    {	  
	      C->Flow_data.Temp_state_vec.Rhof = C->Flow_data.State_vec.Rhof; /*Assign temporary values*/
	      C->Flow_data.Temp_state_vec.RhoU = C->Flow_data.State_vec.RhoU;
	      C->Flow_data.Temp_state_vec.RhoV = C->Flow_data.State_vec.RhoV;
	      C->Flow_data.Temp_state_vec.RhoW = C->Flow_data.State_vec.RhoW;
	      C->Flow_data.Temp_state_vec.RhoE = C->Flow_data.State_vec.RhoE;
	      C->Flow_data.Temp_state_vec.Pres = C->Flow_data.State_vec.Pres;
	      C->Flow_data.Temp_state_vec.T = C->Flow_data.State_vec.T;
     	      
	      if(twoD == 2) {
		area = vol;
#if AXI_OK
		vol = vol*(C->r_centroid);
#else
		vol = vol*(C->centroid[1]); /*Volume per radian*/		
#endif

		/*Extra source term for axisymmetric equations*/
		C -> Flow_data.State_vec.RhoV += ((CFL_data.CFL)*global_dt/(2.0*vol))*(C -> Flow_data.Temp_state_vec.Pres)*area;
	      }
		
	      /*2D or 3D has exact same treatment, except in 2D 'volume' is really area, and 'area' is line length*/
	      C -> Flow_data.State_vec.Rhof -= ((CFL_data.CFL)*global_dt/(2.0*vol))*sum_mss;
	      C -> Flow_data.State_vec.RhoU -= ((CFL_data.CFL)*global_dt/(2.0*vol))*sum_xmnt;
	      C -> Flow_data.State_vec.RhoV -= ((CFL_data.CFL)*global_dt/(2.0*vol))*sum_ymnt;
	      C -> Flow_data.State_vec.RhoW -= ((CFL_data.CFL)*global_dt/(2.0*vol))*sum_zmnt;
	      C -> Flow_data.State_vec.RhoE -= ((CFL_data.CFL)*global_dt/(2.0*vol))*sum_e;

	      if(gas_data.num_species > 1)
		{
		  C->Flow_data.Temp_state_vec.rhof_prod = C->Flow_data.State_vec.rhof_prod;

#if WARHEADBURN
		  if(flowtime <= 0.01e-3) {
#if 0 
 /*Afterburning for helium*/
 C -> Flow_data.State_vec.RhoE += 1.136666667e8*(C->Flow_data.Temp_state_vec.rhof_prod)*((CFL_data.CFL)*global_dt/2.0);
#else
 C -> Flow_data.State_vec.RhoE += 2.04684109763880e+11*(C->Flow_data.Temp_state_vec.rhof_prod)*((CFL_data.CFL)*global_dt/2.0);
#endif

}
#endif

		  C -> Flow_data.State_vec.rhof_prod -= ((CFL_data.CFL)*global_dt/(2.0*vol))*sum_mss_prod;

#if AFTERBURN		  
		  if(ICheat.num_pts != 0) {
		    
		    C->Flow_data.Temp_state_vec.rhof_prod_ub = C->Flow_data.State_vec.rhof_prod_ub;
		    C->Flow_data.State_vec.rhof_prod_ub -= ((CFL_data.CFL)*global_dt/2.0)*sum_mss_prod_ub; 
		    
		    if(C -> Flow_data.State_vec.T >= ICheat.T1) {
		      C -> Flow_data.State_vec.RhoE += (ICheat.e1)*(C->Flow_data.State_vec.rhof_prod_ub)*((CFL_data.CFL)*global_dt/2.0);
		      C->Flow_data.State_vec.rhof_prod_ub -= (ICheat.f1)*((CFL_data.CFL)*global_dt/2.0);

		      /*Order is important - energy released first before products mass fractions are burnt*/
		    }
		    		     
		    if(C->Flow_data.State_vec.rhof_prod_ub <= 0)
		      C->Flow_data.State_vec.rhof_prod_ub = 0; /*Fuel all consumed*/
		  }
#endif
 
		  get_mixture_Pres(&(C -> Flow_data.State_vec)); /*Will get temperature and pressure*/		  
		}
	      else {
		C -> Flow_data.State_vec.Pres = ((gas_data.gam_amb)-1)*
		  ((C->Flow_data.State_vec.RhoE) - 0.5*(SQR(C->Flow_data.State_vec.RhoU) + SQR(C->Flow_data.State_vec.RhoV) +
							SQR(C->Flow_data.State_vec.RhoW))/(C->Flow_data.State_vec.Rhof));
		C -> Flow_data.State_vec.T = (C -> Flow_data.State_vec.Pres)/
		  ((C -> Flow_data.State_vec.Rhof)*((gas_data.Cv_amb)*((gas_data.gam_amb) - 1.0)));
	      }

#if 0
	      printf("C %g %g, sm %g, sx %g, sy %g, se %g, tmp r %g, ru %g, rv %g, re %g, now r %g, ru %g, rv %g, re %g\n",
		     C->centroid[0],C->centroid[1],
		     sum_mss,sum_xmnt,sum_ymnt,sum_e,C->Flow_data.Temp_state_vec.Rhof,C->Flow_data.Temp_state_vec.RhoU,
		     C->Flow_data.Temp_state_vec.RhoV,C->Flow_data.Temp_state_vec.RhoE,C->Flow_data.State_vec.Rhof,
		     C->Flow_data.State_vec.RhoU,C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoE);
#endif

	      if(check_CFL == TRUE)
		{
		  /*In RK step 1 we must see if the CFL no. should be changed*/
		  dr = ABS(((C->Flow_data.State_vec.Rhof) - (C->Flow_data.Temp_state_vec.Rhof))/(C->Flow_data.Temp_state_vec.Rhof));
		  dp = ABS(((C->Flow_data.State_vec.Pres) - (C->Flow_data.Temp_state_vec.Pres))/(C->Flow_data.Temp_state_vec.Pres));

		  /*Get maximum relative changes for this list section*/
		  if(max_dr < dr)
		    max_dr = dr;
		  if(max_dp < dp)
		    max_dp = dp;

		  /*Store flux sums for RK step 1 just in case CFL needs revising*/
		  C -> Flow_data.Flx_sums.sum_mss = sum_mss;
		  C -> Flow_data.Flx_sums.sum_xmnt = sum_xmnt;
		  C -> Flow_data.Flx_sums.sum_ymnt = sum_ymnt;
		  C -> Flow_data.Flx_sums.sum_zmnt = sum_zmnt;
		  C -> Flow_data.Flx_sums.sum_e = sum_e;
		  C -> Flow_data.Flx_sums.sum_mss_prod = sum_mss_prod;
#if AFTERBURN
		  C -> Flow_data.Flx_sums.sum_mss_prod_ub = sum_mss_prod_ub;
#endif
		}	      

#if 0

	      if(C -> Flow_data.State_vec.Pres <= 0) {		
		if(check_CFL == FALSE) {
		  printf("WARNING! trhead %d, P is %e\n",thread_num,C->Flow_data.State_vec.Pres);
		  printf("C %hd, (%e %e %e), type %hd, vol %e\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],C->cell_type,C->cell_volume);
	          printf("thread %d, Temp - r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",thread_num,C->Flow_data.Temp_state_vec.Rhof,
			 C->Flow_data.Temp_state_vec.rhof_prod,C->Flow_data.Temp_state_vec.RhoU,C->Flow_data.Temp_state_vec.RhoV,C->Flow_data.Temp_state_vec.RhoW,C->Flow_data.Temp_state_vec.RhoE,C->Flow_data.Temp_state_vec.Pres);

		  printf("thread %d, r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",thread_num,C->Flow_data.State_vec.Rhof,C->Flow_data.State_vec.rhof_prod,C->Flow_data.State_vec.RhoU,C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE,C->Flow_data.State_vec.Pres);
		  printf("thread %d, CFL now %e, dt %e\n",thread_num,CFL_data.CFL,global_dt);
		  printf("summ %e, sumx %e, sumy %e, sumz %e, sume %e, sumprod %e\n",sum_mss,sum_xmnt,sum_ymnt,sum_zmnt,sum_e,sum_mss_prod);
		  return(ERROR);
		}
	      }
	      else if(C -> Flow_data.State_vec.Pres > 0) {
	      }
	      else if(check_CFL == FALSE) {
		printf("Cell %hd, (%e %e %e)\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
		printf("A %g pressure, RK %d, CFL_cehck %hd!\n",C->Flow_data.State_vec.Pres, RK, check_CFL);
		printf("R %g, ru %g, rv %g, rw %g, re %g\n",C->Flow_data.State_vec.Rhof,C->Flow_data.State_vec.RhoU,
		       C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE);
		printf("Temp - r %g, re %g, p %g\n",C->Flow_data.Temp_state_vec.Rhof,C->Flow_data.Temp_state_vec.RhoE,
		       C->Flow_data.Temp_state_vec.Pres);
		printf("sum_xmnt %g, sum_ymnt %g, sum_zmnt %g\n",sum_xmnt,sum_ymnt,sum_zmnt);
		
		return(ERROR);
	      }
#endif
	      
#if 0
	      if((C -> Flow_data.State_vec.Pres <= 0) && (checking == FALSE)) {
		printf("For cell %hd, (%e %e %e), P %g, RK %d, prev P %g, ic_flag %g, r_p %20.19e\n",C->cell_level,C->centroid[0],
		       C->centroid[1],C->centroid[2],C->Flow_data.State_vec.Pres,RK,C->Flow_data.Temp_state_vec.Pres,C->IC_flag,
		       C->Flow_data.State_vec.rhof_prod);
		printf("sum_mss %g, sum_xmnt %g, sum_ymnt %g, sum_zmnt %g, dt %g, vol %g, ru %g, rv %g, rw %g\n",sum_mss,sum_xmnt,sum_ymnt,sum_zmnt,
		       global_dt,vol,C->Flow_data.State_vec.RhoU, C->Flow_data.State_vec.RhoV, C->Flow_data.State_vec.RhoW);
		printf("re %g, prev_re %g, sum_e %g\n",C->Flow_data.State_vec.RhoE,C->Flow_data.Temp_state_vec.RhoE,sum_e);

#if 1
		if(C->Merge_data!=NULL)
		  printf("Cell's merged\n");
		for(i = 0; i < 6; i++) {
		  printf("Dir %hd - areas ", i);
		  for(j = 0; j < 4; j++)
		    printf("%e ",C -> flux_area[i][j]);
		  printf("\n");
		}
#endif

		for(i = 0; i < 6; i++) {
		  if(((C -> face_neighbours[i][0] == NULL) || (C -> face_neighbours[i][1] == NULL)) && 
		     (C -> Flow_data.Face_fluxes[i][0] != NULL)) {
		    printf("Dir %hd (L %hd) - mf %e, xf %e, yf %e, zf %e, ef %e, mfp %e\n",i,C->face_neighbours[i][0]->cell_level,
			   C->Flow_data.Face_fluxes[i][0]->Mss_flx, C->Flow_data.Face_fluxes[i][0]->X_mntm_flx,
			   C->Flow_data.Face_fluxes[i][0]->Y_mntm_flx,C->Flow_data.Face_fluxes[i][0]->Z_mntm_flx,
			   C->Flow_data.Face_fluxes[i][0]->E_flx,C->Flow_data.Face_fluxes[i][0]->mss_flx_prod);
		  }
		}
		
#if 0
		printf("Ob - mf %e, xf %e, yf %e, zf %e, ef %e, mfp %e\n",C -> Flow_data.Obstructed_flux -> Mss_flx,
		       C -> Flow_data.Obstructed_flux -> X_mntm_flx,C->Flow_data.Obstructed_flux->Y_mntm_flx,
		       C->Flow_data.Obstructed_flux->Z_mntm_flx,C->Flow_data.Obstructed_flux->E_flx,
		       C->Flow_data.Obstructed_flux->mss_flx_prod);
#endif

		return(ERROR);
	      }
#endif

     
#if GET_RID_EPS 
	      if(ABS(C -> Flow_data.State_vec.RhoU) < SMALL_NUM_LIMIT)
		C -> Flow_data.State_vec.RhoU = 0;
	      if(ABS(C -> Flow_data.State_vec.RhoV) < SMALL_NUM_LIMIT)
		C -> Flow_data.State_vec.RhoV = 0;
	      if(ABS(C -> Flow_data.State_vec.RhoW) < SMALL_NUM_LIMIT)
		C -> Flow_data.State_vec.RhoW = 0;
#endif



#if 0
	      if(C -> Flow_data.State_vec.Rhof <= 0) {
		if(check_CFL == FALSE) {
		  printf("RK %d, C->Flow_data.State_vec.Rhof %g\n",RK,C->Flow_data.State_vec.Rhof);
		  printf("Cell's level %hd @ (%g %g %g)\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
		  return(ERROR);
		}
	      } 
	      else if(C -> Flow_data.State_vec.Rhof > 0) {
	      }
	      else {
		if(check_CFL == FALSE) {
		  printf("RK %d,step %d, C->Flow_data.State_vec.Rhof %g, temp Rhof is %g\n",RK,stepnow,C->Flow_data.State_vec.Rhof,
			 C->Flow_data.Temp_state_vec.Rhof);
		  printf("C L %hd, (%e %e %e)\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
		  printf("sum_mss is %g, global_dt %g, CFL %g, rc %g, yc %g\n",sum_mss,global_dt,CFL_data.CFL,C->r_centroid,
			 C->centroid[1]);
		  if(C->Merge_data != NULL)
		    printf("C's merged\n");
		  return(ERROR);
		}
	      }

	      if((C -> Flow_data.State_vec.RhoE <= 0) && (check_CFL == FALSE)) {
		printf("C->Flow_data.State_vec.RhoE %g\n",C->Flow_data.State_vec.RhoE);
		return(ERROR);
	      } 
	      else if(C -> Flow_data.State_vec.RhoE > 0) {
	      }
	      else {
		if(check_CFL == FALSE) {
		  printf("C->Flow_data.State_vec.RhoE %g\n",C->Flow_data.State_vec.RhoE);
		  return(ERROR);
		}
	      }
#endif

	      if(C -> Merge_data != NULL) 
		{ /*Merged cells also take on normal cell's values - we need only to
		    store the temporary state vectors within the normal cell as all cells within this cluster have this
		    same value; updating small cells simply involves looking at the normal cell's values*/
		  
		  merge_node = C -> Merge_data -> Linked_cells;
		  
		  while(merge_node != NULL)
		    {
		      S = merge_node -> cell_loc;
		      
		      if(S -> is_small_cell == TRUE)
			{
			  S->Flow_data.State_vec.Rhof = C->Flow_data.State_vec.Rhof; 
			  S->Flow_data.State_vec.RhoU = C->Flow_data.State_vec.RhoU;
			  S->Flow_data.State_vec.RhoV = C->Flow_data.State_vec.RhoV;
			  S->Flow_data.State_vec.RhoW = C->Flow_data.State_vec.RhoW;
			  S->Flow_data.State_vec.RhoE = C->Flow_data.State_vec.RhoE;
			  S->Flow_data.State_vec.Pres = C->Flow_data.State_vec.Pres;
			  S->Flow_data.State_vec.T = C->Flow_data.State_vec.T;

			  if(gas_data.num_species > 1) {
			    S->Flow_data.State_vec.rhof_prod = C->Flow_data.State_vec.rhof_prod;
#if AFTERBURN
			    S->Flow_data.State_vec.rhof_prod_ub = C->Flow_data.State_vec.rhof_prod_ub;
#endif
			  }			    
			} 
		      
		      merge_node = merge_node -> next;
		    }
		}		    
	    }
	  else if(stage == 2) /*2nd stage*/
	    {	  
	      if(twoD == 2) {
		temp_pres = C -> Flow_data.State_vec.Pres;
		area = vol;

#if AXI_OK
		vol = vol*(C->r_centroid);
#else
		vol = vol*(C->centroid[1]);		
#endif
	      }

	      C->Flow_data.State_vec.Rhof = (C->Flow_data.Temp_state_vec.Rhof)-((CFL_data.CFL)*global_dt/vol)*sum_mss;
	      C->Flow_data.State_vec.RhoU = (C->Flow_data.Temp_state_vec.RhoU)-((CFL_data.CFL)*global_dt/vol)*sum_xmnt;
	      C->Flow_data.State_vec.RhoV = (C->Flow_data.Temp_state_vec.RhoV)-((CFL_data.CFL)*global_dt/vol)*sum_ymnt;
	      C->Flow_data.State_vec.RhoW = (C->Flow_data.Temp_state_vec.RhoW)-((CFL_data.CFL)*global_dt/vol)*sum_zmnt;
	      C->Flow_data.State_vec.RhoE = (C->Flow_data.Temp_state_vec.RhoE)-((CFL_data.CFL)*global_dt/vol)*sum_e;

	      if(twoD == 2)
		C -> Flow_data.State_vec.RhoV += ((CFL_data.CFL)*global_dt/vol)*temp_pres*area;
	      
	      if(gas_data.num_species > 1)
		{
#if WARHEADBURN
                  if(flowtime <= 0.01e-3) {
#if 0
 C -> Flow_data.State_vec.RhoE += 1.136666667e8*(C->Flow_data.State_vec.rhof_prod)*((CFL_data.CFL)*global_dt);
#else
 C -> Flow_data.State_vec.RhoE += 2.04684109763880e+11*(C->Flow_data.State_vec.rhof_prod)*((CFL_data.CFL)*global_dt); 
#endif
}
#endif

		  C->Flow_data.State_vec.rhof_prod = (C->Flow_data.Temp_state_vec.rhof_prod)-((CFL_data.CFL)*global_dt/vol)*sum_mss_prod;

#if AFTERBURN
		  if(ICheat.num_pts != 0) {
		    C->Flow_data.State_vec.rhof_prod_ub = (C->Flow_data.Temp_state_vec.rhof_prod_ub) - 
		      sum_mss_prod_ub*((CFL_data.CFL)*global_dt);

		    if(C -> Flow_data.Temp_state_vec.T >= ICheat.T2) {
		      C -> Flow_data.State_vec.RhoE += (ICheat.e2)*(C->Flow_data.State_vec.rhof_prod_ub)*((CFL_data.CFL)*global_dt);
		      C->Flow_data.State_vec.rhof_prod_ub -= (ICheat.f2)*((CFL_data.CFL)*global_dt);
		    }

		    if(C->Flow_data.State_vec.rhof_prod_ub <= 0)
		      C->Flow_data.State_vec.rhof_prod_ub = 0;
		  }
#endif

		  get_mixture_Pres(&(C -> Flow_data.State_vec)); /*Will get temperature and pressure*/
		}

	      else {
		C -> Flow_data.State_vec.Pres = ((gas_data.gam_amb)-1.0)*
		  ((C->Flow_data.State_vec.RhoE) - 0.5*(SQR(C->Flow_data.State_vec.RhoU) + SQR(C->Flow_data.State_vec.RhoV) +
							SQR(C->Flow_data.State_vec.RhoW))/(C->Flow_data.State_vec.Rhof));
		C -> Flow_data.State_vec.T = (C -> Flow_data.State_vec.Pres)/
		  ((C -> Flow_data.State_vec.Rhof)*((gas_data.Cv_amb)*((gas_data.gam_amb) - 1.0)));
	      }

#if GET_RID_EPS 
	      if(ABS(C -> Flow_data.State_vec.RhoU) < SMALL_NUM_LIMIT)
		C -> Flow_data.State_vec.RhoU = 0;
	      if(ABS(C -> Flow_data.State_vec.RhoV) < SMALL_NUM_LIMIT)
		C -> Flow_data.State_vec.RhoV = 0;
	      if(ABS(C -> Flow_data.State_vec.RhoW) < SMALL_NUM_LIMIT)
		C -> Flow_data.State_vec.RhoW = 0;
#endif

#if 0

	      if(C -> Flow_data.State_vec.Pres <= 0) {		
		if(check_CFL == FALSE) {
		  printf("WARNING! trhead %d, P is %e\n",thread_num,C->Flow_data.State_vec.Pres);
		  printf("C %hd, (%e %e %e), type %hd, vol %e\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],C->cell_type,C->cell_volume);
	          printf("thread %d, Temp - r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",thread_num,C->Flow_data.Temp_state_vec.Rhof,
			 C->Flow_data.Temp_state_vec.rhof_prod,C->Flow_data.Temp_state_vec.RhoU,C->Flow_data.Temp_state_vec.RhoV,C->Flow_data.Temp_state_vec.RhoW,C->Flow_data.Temp_state_vec.RhoE,C->Flow_data.Temp_state_vec.Pres);
		  printf("RK0.5 pres is %g\n",temp_pres);

		  printf("thread %d, r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",thread_num,C->Flow_data.State_vec.Rhof,C->Flow_data.State_vec.rhof_prod,C->Flow_data.State_vec.RhoU,C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE,C->Flow_data.State_vec.Pres);
		  printf("thread %d, CFL now %e, dt %e\n",thread_num,CFL_data.CFL,global_dt);
		  printf("summ %e, sumx %e, sumy %e, sumz %e, sume %e, sumprod %e\n",sum_mss,sum_xmnt,sum_ymnt,sum_zmnt,sum_e,sum_mss_prod);
		  return(ERROR);
		}
	      }
	      else if(C -> Flow_data.State_vec.Pres > 0) {
	      }
	      else if(check_CFL == FALSE) {
		printf("Cell %hd, (%e %e %e)\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
		printf("A %g pressure, RK %d, CFL_cehck %hd!\n",C->Flow_data.State_vec.Pres, RK, check_CFL);
		printf("R %g, ru %g, rv %g, rw %g, re %g\n",C->Flow_data.State_vec.Rhof,C->Flow_data.State_vec.RhoU,
		       C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE);
		printf("Temp - r %g, re %g, p %g\n",C->Flow_data.Temp_state_vec.Rhof,C->Flow_data.Temp_state_vec.RhoE,
		       C->Flow_data.Temp_state_vec.Pres);
		printf("sum_xmnt %g, sum_ymnt %g, sum_zmnt %g\n",sum_xmnt,sum_ymnt,sum_zmnt);
		
		return(ERROR);
	      }


	      if(C -> Flow_data.State_vec.Rhof <= 0) {
		printf("RK %d,C->Flow_data.State_vec.Rhof %g\n",RK,C->Flow_data.State_vec.Rhof);
		return(ERROR);
	      } 
	      else if(C -> Flow_data.State_vec.Rhof > 0) {
	      }
	      else {
		printf("RK %d,C->Flow_data.State_vec.Rhof %g\n",RK,C->Flow_data.State_vec.Rhof);
		return(ERROR);
	      }

	      if(C -> Flow_data.State_vec.RhoE <= 0) {
		printf("C->Flow_data.State_vec.RhoE %g\n",C->Flow_data.State_vec.RhoE);
		return(ERROR);
	      } 
	      else if(C -> Flow_data.State_vec.RhoE > 0) {
	      }
	      else {
		printf("C->Flow_data.State_vec.RhoE %g\n",C->Flow_data.State_vec.RhoE);
		return(ERROR);
	      }

	      if(C -> Flow_data.State_vec.Pres <= 0) {
		printf("RK %d, C->Flow_data.State_vec.Pres %g\n",RK, C->Flow_data.State_vec.Pres);
		printf("C's cell type is %hd, level %hd, @ (%g %g %g)\n",C->cell_type,C->cell_level,C->centroid[0],C->centroid[1],
		       C->centroid[2]);
		return(ERROR);
	      } 
	      else if(C -> Flow_data.State_vec.Pres > 0) {
	      }
	      else {
		printf("C->Flow_data.State_vec.Pres %g\n",C->Flow_data.State_vec.Pres);
		return(ERROR);
	      }
#endif

#if 0
	      sc_count = 0;
#endif
	      if(C -> Merge_data != NULL)
		{
		  merge_node = C -> Merge_data -> Linked_cells;
		  
		  while(merge_node != NULL)
		    {
		      S = merge_node -> cell_loc;
#if 0
		      sc_count++;
#endif
		      if(S -> is_small_cell == TRUE)
			{
			  S->Flow_data.State_vec.Rhof = C->Flow_data.State_vec.Rhof; 
			  S->Flow_data.State_vec.RhoU = C->Flow_data.State_vec.RhoU;
			  S->Flow_data.State_vec.RhoV = C->Flow_data.State_vec.RhoV;
			  S->Flow_data.State_vec.RhoW = C->Flow_data.State_vec.RhoW;
			  S->Flow_data.State_vec.RhoE = C->Flow_data.State_vec.RhoE;
			  S->Flow_data.State_vec.Pres = C->Flow_data.State_vec.Pres;
			  S->Flow_data.State_vec.T = C->Flow_data.State_vec.T;

			  if(gas_data.num_species > 1) {
			    S->Flow_data.State_vec.rhof_prod = C->Flow_data.State_vec.rhof_prod;
#if AFTERBURN
			    S->Flow_data.State_vec.rhof_prod_ub = C->Flow_data.State_vec.rhof_prod_ub;
#endif
			  }			    
			}
		      
		      merge_node = merge_node -> next;
		    }
#if 0
		  if(sc_count == 2)
		    printf("sc_count is 2, thread %d\n",omp_get_thread_num());
		  else if(sc_count >= 3)
		    printf("More than 3, thread %d\n",omp_get_thread_num());
#endif
		}	  
	    }  

	  if(twoD != 0)
	    C -> Flow_data.State_vec.RhoW = 0; /*Extra precautionary measure*/
	    
	}
      
      if(node == Tail) /*Only look at list region between Head and Tail*/
	break;
      else node = node -> next;
    }  

  if(check_CFL == TRUE) /*Need only to revise CFL once for RK stage 1, since it'll be revised significantly downward*/
    {
      CFL_data.CFL_min[thread_num] = MIN(((CFL_data.CFL_cutback)/MAX(max_dr, max_dp)), CFL_data.CFL_max);
      if(CFL_data.CFL_min[thread_num] < CFL_data.CFL_max) 
	CFL_data.CFL_ok = FALSE; /*CFL no. must be revised; previous CFL is maximum CFL*/
     
      #pragma omp barrier
      /*Must unsure all threads reach this point so CFL_ok will be known by all*/

      if(CFL_data.CFL_ok == FALSE)
	{ /*We must reset the values of the state vectors (only necessary if we need to undergo time integration again)*/

	  node = L;
	  while(node != NULL) 
	    {
	      C = node -> cell_loc;

	      if(
#if DETON
		 (C -> un_det != TRUE) &&
#endif
		 (C -> is_small_cell == FALSE))
		{
		  C -> Flow_data.State_vec.Rhof = C -> Flow_data.Temp_state_vec.Rhof;
		  C -> Flow_data.State_vec.RhoU = C -> Flow_data.Temp_state_vec.RhoU;
		  C -> Flow_data.State_vec.RhoV = C -> Flow_data.Temp_state_vec.RhoV;
		  C -> Flow_data.State_vec.RhoW = C -> Flow_data.Temp_state_vec.RhoW;
		  C -> Flow_data.State_vec.RhoE = C -> Flow_data.Temp_state_vec.RhoE;
		  C -> Flow_data.State_vec.Pres = C -> Flow_data.Temp_state_vec.Pres;
		  C -> Flow_data.State_vec.T = C -> Flow_data.Temp_state_vec.T;
		  C -> Flow_data.State_vec.rhof_prod = C -> Flow_data.Temp_state_vec.rhof_prod;
#if AFTERBURN
		  C->Flow_data.State_vec.rhof_prod_ub = C -> Flow_data.Temp_state_vec.rhof_prod_ub;
#endif

		  if(C -> Merge_data != NULL)
		    {
		      merge_node = C -> Merge_data -> Linked_cells;

		      while(merge_node != NULL)
			{
			  S = merge_node -> cell_loc;

			  if(S -> is_small_cell == TRUE)
			    {
			      S -> Flow_data.State_vec.Rhof = C -> Flow_data.Temp_state_vec.Rhof;
			      S -> Flow_data.State_vec.RhoU = C -> Flow_data.Temp_state_vec.RhoU;
			      S -> Flow_data.State_vec.RhoV = C -> Flow_data.Temp_state_vec.RhoV;
			      S -> Flow_data.State_vec.RhoW = C -> Flow_data.Temp_state_vec.RhoW;
			      S -> Flow_data.State_vec.RhoE = C -> Flow_data.Temp_state_vec.RhoE;
			      S -> Flow_data.State_vec.Pres = C -> Flow_data.Temp_state_vec.Pres;
			      S -> Flow_data.State_vec.T = C -> Flow_data.Temp_state_vec.T;
			      S -> Flow_data.State_vec.rhof_prod = C -> Flow_data.Temp_state_vec.rhof_prod;
#if AFTERBURN
			      S->Flow_data.State_vec.rhof_prod_ub = C -> Flow_data.Temp_state_vec.rhof_prod_ub;
#endif
			    }

			  merge_node = merge_node -> next; 
			}
		    }
		}

	      if(node == Tail) 
		break;
	      else node = node -> next;
	    }	  
	}
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Traverse through list of cells and exchange fluxes (i.e. calculate what the flux vectors are for each face/quadrant).
   These fluxes values will then be used later (summed over all areas) for time integration.*/

short int all_cells_exchange_fluxes(List_leaf L, List_leaf Tail, double time, short int stage)
{
  List_leaf node = L;
  Cart_cell C;
  State_vector boundary_state, Lstate, Rstate;
  short int i, j, k, opp_dir, reconstructed_yet; 
  short int dir[3] = {EST, NTH, UPR}; /*'Positive' directions*/
  double norm[3];
  char recon_full;

#if DETON
  double flux_area[4], flux_area_opp[4];
#endif
  opp_dir = 0;
  while(node != NULL) /*Traverse through list*/
    {
      C = node -> cell_loc;

#if DETON
      if(C -> un_det != TRUE)
	{
#endif
	  if(wall_rep == 1) /*As only 1 'obstructed area', can compute flux here*/
	    {
	      for(j = 0; j < 3; j++) 
		{
		  i = 2*j;
		  opp_dir = i+1;
	  
#if DETON
		  if(Deton.num_deton_pts > 0) {
		    for(k = 0; k < 4; k++) { /*First get the 'proper' flux areas*/
		      flux_area[k] = C->flux_area[i][k];
		      flux_area_opp[k] = C->flux_area[opp_dir][k];
		    }

		    /*Now change the flux areas accordingly*/
		    if(C -> face_neighbours[i][0] != NULL) {
		      if((C -> face_neighbours[i][1] == NULL)) {
			if(C -> face_neighbours[i][0] -> un_det == TRUE) {		      
			  for(k = 0; k < 4; k++)
			    flux_area[k] = 0; /*Undetonated cell is like a solid cell*/
			}
		      }
		      else {
			for(k = 0; k < 4; k++) {
			  if(C -> face_neighbours[i][k] -> un_det == TRUE)
			    flux_area[k] = 0;
			}
		      }
		    }

		    if(C -> face_neighbours[opp_dir][0] != NULL) {
		      if((C -> face_neighbours[opp_dir][1] == NULL)) {
			if(C -> face_neighbours[opp_dir][0] -> un_det == TRUE) {		      
			  for(k = 0; k < 4; k++)
			    flux_area_opp[k] = 0; 
			}
		      }
		      else {
			for(k = 0; k < 4; k++) {
			  if(C -> face_neighbours[opp_dir][k] -> un_det == TRUE)
			    flux_area_opp[k] = 0;
			}
		      }
		    }

		    C->wall_norm[j] = (flux_area_opp[0]+flux_area_opp[1]+flux_area_opp[2]+flux_area_opp[3]) - 
		      (flux_area[0]+flux_area[1]+flux_area[2]+flux_area[3]);
		  }
		  else {
		    C->wall_norm[j] = (C->flux_area[opp_dir][0]+C->flux_area[opp_dir][1]+C->flux_area[opp_dir][2]+C->flux_area[opp_dir][3])
		      - (C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]);
		  } 
#else
		  C -> wall_norm[j] = (C->flux_area[opp_dir][0]+C->flux_area[opp_dir][1]+C->flux_area[opp_dir][2]+C->flux_area[opp_dir][3])
		    - (C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]);
#endif
		}
      
	      if((ABS(C->wall_norm[0]) > AREA_DIFF_TOL) || (ABS(C->wall_norm[1]) > AREA_DIFF_TOL) || 
		 (ABS(C->wall_norm[2]) > AREA_DIFF_TOL)) /*An obstruction exists*/
		{
		  C -> wall_area = SQRT(SQR(C->wall_norm[0])+SQR(C->wall_norm[1])+SQR(C->wall_norm[2])); /*Area is magnitude of normal*/

		  for(j = 0; j < 3; j++)
		    C -> wall_norm[j] = ((C -> wall_norm[j])/(C -> wall_area)); /*Normalize*/

		  if(C -> Flow_data.Obstructed_flux == NULL)
		    {
		      C -> Flow_data.Obstructed_flux = malloc(sizeof(struct flux_vector));
		      if(C -> Flow_data.Obstructed_flux == NULL)
			return(ERROR);
		    }

		  if((twoD == 2) && (C->cell_type == FLUID) && (((ABS(C->wall_norm[0]) < 1) && (ABS(C->wall_norm[0]) > AREA_DIFF_TOL)) && 
								((ABS(C->wall_norm[1]) < 1) && (ABS(C->wall_norm[1]) > AREA_DIFF_TOL)))) {
		    /*Fluid cells with wall norms that aren't axis unit vectors i.e. grid-aligned cells, are better
		      treated as staircases, but do this only for axisymmetric (because of source term)*/
		  }
		  else {
		    boundary_state = set_boundary_state(&(C->Flow_data.State_vec), NULL, C->wall_norm, IRRELEVANT, WALL,
							IRRELEVANT, IRRELEVANT, IRRELEVANT, NULL,0);
		    if(flux_type == ADAPTIVE) {
		      if(C -> shocked_cell == 'y')
			compute_fluxes(&(C->Flow_data.State_vec), &boundary_state, C->wall_norm, C->Flow_data.Obstructed_flux, EFM);
		      else compute_fluxes(&(C->Flow_data.State_vec), &boundary_state, C->wall_norm, C->Flow_data.Obstructed_flux, AUSMDV);
		    }
		    else compute_fluxes(&(C->Flow_data.State_vec), &boundary_state, C->wall_norm, C->Flow_data.Obstructed_flux, flux_type);
		  } 
		}
	      else 
		{
		  C -> wall_area = 0;
		  C -> wall_norm[0] = 0;
		  C -> wall_norm[1] = 0;
		  C -> wall_norm[2] = 0;
		}
	    }
            
	  /*For boundary cells use full (unlimited) reconstruction (EFM will also be used),
	   as non-reflecting BCs (if used) do better*/

	  recon_full = 'n';
	  if(want_to_adapt == FALSE) {
	    for(i = 0; i < NUM_FACES; i++) {
	      if(C -> face_neighbours[i][0] == NULL) {
		recon_full = 'y';
		break;
	      }
	    }
	  }

	  for(i = 0; i < NUM_FACES; i++) /*First try and compute border fluxes*/
	    {
	      norm[0] = 0; norm[1] = 0; norm[2] = 0;

	      switch(i)
		{
		case EST:
		  norm[0] = 1; break;
		case WST:
		  norm[0] = -1; break;
		case NTH:
		  norm[1] = 1; break;
		case STH:
		  norm[1] = -1; break;
		case UPR:
		  norm[2] = 1; break;
		case LWR:
		  norm[2] = -1; break;	      
		}

	      if((C -> face_neighbours[i][0] == NULL) && (C -> Flow_data.Face_fluxes[i][0] != NULL)) 
		{
		  Lstate = reconstruct_flow(C, i, IRRELEVANT, IRRELEVANT, IRRELEVANT, EOF, NULL, NULL, recon_full);
		  /*We basically want to extrapolate i.e. use limiter value of 1.0*/

#if SUPERSONIC_VORTEX
		  Lstate.r = SQRT(SQR(C->centroid[0]) + SQR(C->centroid[1]));
#endif

		  /*Transient values no need to set the boundary state for individual cells?*/
		  if(many_limit == 'n')
		    boundary_state = set_boundary_state(&Lstate, &(C->Flow_data.State_vec), NULL, i, BORDER, stage, time, 
							0.5*CALC_CELL_EDGE_LENGTH(C), &(C -> Flow_data.Grads), C->Flow_data.Lim.lim);
		  else {
		    if((i == EST) || (i == WST)) 
		      boundary_state = set_boundary_state(&Lstate, &(C->Flow_data.State_vec), NULL, i, BORDER, stage, time, 
							  0.5*CALC_CELL_EDGE_LENGTH(C), &(C -> Flow_data.Grads), C->Flow_data.Lim.u_lim);
		    else if((i == NTH) || (i == STH)) 
		      boundary_state = set_boundary_state(&Lstate, &(C->Flow_data.State_vec), NULL, i, BORDER, stage, time, 
							  0.5*CALC_CELL_EDGE_LENGTH(C), &(C -> Flow_data.Grads), C->Flow_data.Lim.v_lim);
		    else if((i == UPR) || (i == LWR)) 
		      boundary_state = set_boundary_state(&Lstate, &(C->Flow_data.State_vec), NULL, i, BORDER, stage, time, 
							  0.5*CALC_CELL_EDGE_LENGTH(C), &(C -> Flow_data.Grads), C->Flow_data.Lim.w_lim);
		  }

#if 0
		  if(((boundary_state.RhoE) - 
		      0.5*(SQR(boundary_state.RhoU) + SQR(boundary_state.RhoV) + SQR(boundary_state.RhoW))/(boundary_state.Rhof)) < 0) {
		    printf("It's boundary here\n");
		    return(ERROR);
		  }
#endif

#if 0
		  if(Lstate.Rhof <= 0) {
		    printf("Lstate.Rhof <= 0, it's %g, boundary dir %hd\n",Lstate.Rhof, i);
		    printf("C's cell type %hd, vol %g, cen %9.8e %9.8e %9.8e\n",C->cell_type,C->cell_volume,C->centroid[0],C->centroid[1],
			   C->centroid[2]);
		    printf("Flux areas - \n");
		    for(i = 0; i < 6; i++)
		      printf("Dir %hd, %e %e %e %e\n",i,C->flux_area[i][0],C->flux_area[i][1],C->flux_area[i][2],C->flux_area[i][3]);
		    printf("Other Lstate values - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",Lstate.Rhof,Lstate.rhof_prod,
			   Lstate.RhoU,Lstate.RhoV,Lstate.RhoW,Lstate.RhoE,Lstate.Pres);
		    printf("C's values - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->Flow_data.State_vec.Rhof,C->Flow_data.State_vec.rhof_prod,C->Flow_data.State_vec.RhoU,C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE,C->Flow_data.State_vec.Pres);
		    printf("C's r grad - %e %e %e\n",C->Flow_data.Grads.grad_rhof[0],C->Flow_data.Grads.grad_rhof[1],C->Flow_data.Grads.grad_rhof[2]);
		    printf("C's limit %e\n",C->Flow_data.Lim.lim);
		    printf("Neighbours' values - \n");
		    for(i = 0; i < 6; i++) {
		      if((C->face_neighbours[i][0] != NULL) && (C->face_neighbours[i][1] == NULL)) {
			printf("1 neighb dir %hd (L %hd, type %hd) - r %g\n",i,C->face_neighbours[i][0]->cell_level,
			       C->face_neighbours[i][0]->cell_type,C->face_neighbours[i][0]->Flow_data.State_vec.Rhof);
		      }
		      else if(C->face_neighbours[i][1] != NULL) {
			for(j = 0; j < 4; j++) {
			  printf("4 neighb dir %hd quad %hd (type %hd) - r %g\n",i,j,C->face_neighbours[i][j]->cell_type,
				 C->face_neighbours[i][j]->Flow_data.State_vec.Rhof);
			}
		      }
		    }
		    return(ERROR);
		  }
		  else if(Lstate.Rhof > 0) {
		  }
		  else {
		    printf("Lstate Rhof bad - %g, dir %hd, RK %d\n",Lstate.Rhof,i,RK);
		    return(ERROR);
		  }
#endif
#if 0
		  if(boundary_state.Rhof <= 0) {
		    printf("boundary_state.Rhof <= 0, it's %g, boundary\n",boundary_state.Rhof);
		    return(ERROR);
		  }
		  else if(boundary_state.Rhof > 0) {
		  }
		  else {
		    printf("boundary Rhof bad boundary - %g, dir %hd, RK %d\n",boundary_state.Rhof, i,RK);
		    return(ERROR);
		  }

		  if(boundary_state.Pres <= 0) {
		    printf("boundary_state.Pres <= 0, it's %g, boundary\n",boundary_state.Pres);
		    printf("Cell in question is level %hd @ (%g %g %g)\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
		    return(ERROR);
		  }
		  else if(boundary_state.Pres > 0) {
		  }
		  else {
		    printf("boundary Pres bad boundary - %g, dir %hd, RK %d\n",boundary_state.Pres, i,RK);
		    return(ERROR);
		  }
#endif

		  /*Always use EFM at boundary as reflections are minimized*/
		  compute_fluxes(&Lstate, &boundary_state, norm, C -> Flow_data.Face_fluxes[i][0], EFM);
		}
	    }
      
	  for(i = 0; i < 3; i++) /*Now compute cell to cell fluxes only in 'positive' directions*/
	    {
	      switch(dir[i])
		{
		case EST:
		  norm[0] = 1; norm[1] = 0; norm[2] = 0;
		  opp_dir = WST; break;
		case NTH:
		  norm[0] = 0; norm[1] = 1; norm[2] = 0;
		  opp_dir = STH; break;
		case UPR:
		  norm[0] = 0; norm[1] = 0; norm[2] = 1;
		  opp_dir = LWR; break;
		}
	  
	      if((C -> face_neighbours[dir[i]][0] != NULL) && (C -> face_neighbours[dir[i]][1] == NULL) && 
#if DETON
		 (C -> face_neighbours[dir[i]][0] -> un_det != TRUE) &&
#endif
		 (C -> Flow_data.Face_fluxes[dir[i]][0] != NULL))
		{ /*Only 1 neighbour and can flux through face*/

		  Lstate = reconstruct_flow(C, dir[i], IRRELEVANT, IRRELEVANT, IRRELEVANT, EOF, NULL, NULL, recon_full);

		  if((C -> face_neighbours[dir[i]][0] -> cell_level) < (C -> cell_level)) /*1 level coarser neighbour*/
		    {
		      switch(dir[i])
			{
			case EST:
			  if(C -> child_num <= 3)
			    Rstate = reconstruct_flow(C -> face_neighbours[dir[i]][0], opp_dir, (C->child_num)-2, 
						      IRRELEVANT, IRRELEVANT, EOF, NULL, NULL, 'n'); /*Quadrants of neighbour on 
												       opposite face*/
			  else Rstate = reconstruct_flow(C -> face_neighbours[dir[i]][0], opp_dir, (C->child_num)-4, 
							 IRRELEVANT, IRRELEVANT, EOF, NULL, NULL, 'n');
			  break;
			case NTH:
			  Rstate = reconstruct_flow(C -> face_neighbours[dir[i]][0], opp_dir, ((C->child_num)-1)/2, 
						    IRRELEVANT, IRRELEVANT, EOF, NULL, NULL, 'n');
			  break;
			case UPR:
			  Rstate = reconstruct_flow(C -> face_neighbours[dir[i]][0], opp_dir, (C->child_num)-4, 
						    IRRELEVANT, IRRELEVANT, EOF, NULL, NULL, 'n');
			  break;
			}
		    }
		  else Rstate = reconstruct_flow(C -> face_neighbours[dir[i]][0], opp_dir, 
						 IRRELEVANT, IRRELEVANT, IRRELEVANT, EOF, NULL, NULL, 'n'); /*Same level neighbour*/


#if 0
		  if(((Rstate.RhoE) - 0.5*(SQR(Rstate.RhoU) + SQR(Rstate.RhoV) + SQR(Rstate.RhoW))/(Rstate.Rhof)) < 0) {
		    printf("Will get a -ve energy here of %g, dir %hd\n", ((Rstate.RhoE) - 0.5*(SQR(Rstate.RhoU) + SQR(Rstate.RhoV) + SQR(Rstate.RhoW))/(Rstate.Rhof))/(Rstate.Rhof), dir[i]);
		    printf("C %hd (%9.8e %9.8e %9.8e), child_num %hd\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],C->child_num);
		    printf("type %hd, vol %e\n",C->cell_type,C->cell_volume);
		    printf("C -> r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",C->Flow_data.State_vec.Rhof,C->Flow_data.State_vec.rhof_prod,
			   C->Flow_data.State_vec.RhoU,C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE,
			   C->Flow_data.State_vec.Pres);
		    printf("Lstate -> r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",Lstate.Rhof,Lstate.rhof_prod,Lstate.RhoU,
			   Lstate.RhoV,Lstate.RhoW,Lstate.RhoE,Lstate.Pres);
		
		    printf("Neighb %hd (%9.8e %9.8e %9.8e),type %hd, vol %e\n",C->face_neighbours[dir[i]][0]->cell_level,C->face_neighbours[dir[i]][0]->centroid[0],
			   C->face_neighbours[dir[i]][0]->centroid[1],C->face_neighbours[dir[i]][0]->centroid[2],C->face_neighbours[dir[i]][0]->cell_type,C->face_neighbours[dir[i]][0]->cell_volume);
		    
		    printf("Neighb -> r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Rhof,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.rhof_prod,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoU,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoV,C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoW,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoE,C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Pres);
		    printf("Neighb -> temp r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.Rhof,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.rhof_prod,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoU,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoV,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoW,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoE,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.Pres);

		    printf("Rstate -> r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e\n",Rstate.Rhof,Rstate.rhof_prod,Rstate.RhoU,
			   Rstate.RhoV,Rstate.RhoW,Rstate.RhoE,Rstate.Pres);

		    return(ERROR);
		  }
#endif

#if 0
		  if(Lstate.Rhof <= 0) { 
		    printf("Lstate.Rhof <= 0, it's %g, 1 neighb dir %hd\n",Lstate.Rhof, i);
		    return(ERROR);
		  }
		  else if(Lstate.Rhof > 0) {
		  }
		  else {
		    printf("Lstate Rhof bad - %g, dir %hd, RK %d\n",Lstate.Rhof,dir[i],RK);
		    return(ERROR);
		  }
#endif
#if 0
		  if(Rstate.Rhof <= 0) {
		    printf("Rstate.Rhof <= 0, it's %g, 1 neighb dir %hd, RK %d\n",Rstate.Rhof, i, RK);
		    return(ERROR);
		  }
		  else if(Rstate.Rhof > 0) {
		  }
		  else {
		    printf("Rstate Rhof bad (L %hd (%e %e %e), type %hd) - %g, dir %hd, RK %d\n",C->cell_level,
			   C->centroid[0],C->centroid[1],C->centroid[2],C->cell_type,Rstate.Rhof, dir[i],RK);
		    printf("But neigbh rho is %g\n",C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Rhof);
		    printf("Neighb L %hd, type %hd\n",C->face_neighbours[dir[i]][0]->cell_level,C->face_neighbours[dir[i]][0]->cell_type);
		    return(ERROR);
		  }

		  if(Rstate.Pres <= 0) { 
		    printf("RK %d, Rstate.Pres <= 0, it's %g, 1 neighb dir %hd\n",RK, Rstate.Pres, dir[i]);
		    printf("Cell %hd (%e %e %e), type %hd, vol %g\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],
			   C->cell_type,C->cell_volume);
		    printf("Cell - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->Flow_data.State_vec.Rhof,
			   C->Flow_data.State_vec.rhof_prod,C->Flow_data.State_vec.RhoU,C->Flow_data.State_vec.RhoV,
			   C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE,C->Flow_data.State_vec.Pres);
		    printf("Cell prev - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->Flow_data.Temp_state_vec.Rhof,
			   C->Flow_data.Temp_state_vec.rhof_prod,C->Flow_data.Temp_state_vec.RhoU,C->Flow_data.Temp_state_vec.RhoV,
			   C->Flow_data.Temp_state_vec.RhoW,C->Flow_data.Temp_state_vec.RhoE,C->Flow_data.Temp_state_vec.Pres);
		    printf("Neigb - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Rhof,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.rhof_prod,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoU,C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoV,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoW,C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoE,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Pres);
		    printf("Neigb prev - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.Rhof,
			   C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.rhof_prod,
			   C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoU,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoV,
			   C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoW,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoE,
			   C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.Pres);
		    printf("Neighbour L %hd (%e %e %e), type %hd, vol %g, Pres %g\n",C->face_neighbours[dir[i]][0]->cell_level,
			   C->face_neighbours[dir[i]][0]->centroid[0],C->face_neighbours[dir[i]][0]->centroid[1],
			   C->face_neighbours[dir[i]][0]->centroid[2],
			   C->face_neighbours[dir[i]][0]->cell_type,C->face_neighbours[dir[i]][0]->cell_volume,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Pres);
		    printf("%hd flx areas - %g %g %g %g\n",dir[i],C->flux_area[dir[i]][0],C->flux_area[dir[i]][1],C->flux_area[dir[i]][2],
			   C->flux_area[dir[i]][3]);
		    printf("opp_dir %hd flx areas - %g %g %g %g\n",opp_dir,C->face_neighbours[dir[i]][0]->flux_area[opp_dir][0],
			   C->face_neighbours[dir[i]][0]->flux_area[opp_dir][1],C->face_neighbours[dir[i]][0]->flux_area[opp_dir][2],
			   C->face_neighbours[dir[i]][0]->flux_area[opp_dir][3]);
		    
		    printf("Rstate - r %g, ru %g, rv %g, rw %g, re %g, rp %g\n",Rstate.Rhof,Rstate.RhoU,Rstate.RhoV,Rstate.RhoW,
			   Rstate.RhoE, Rstate.rhof_prod);
		    if(C->Merge_data != NULL) {
		      printf("This cell's merged\n");
		      if(C->Merge_data == C->face_neighbours[dir[i]][0]->Merge_data)
			printf("Merged with neighb\n");
		    }

		    if(C->face_neighbours[dir[i]][0]->Merge_data != NULL)
		      printf("Neighbour's also merged\n");
		    if(C->is_small_cell == TRUE)
		      printf("C's a small cell\n");
		    if(C->face_neighbours[dir[i]][0]->is_small_cell == TRUE)
		      printf("Neighb's a small cell\n");

		    return(ERROR);
		  }
		  else if(Rstate.Pres > 0) {
		  }
		  else {
		    printf("Rstate Pres bad - %g, dir %hd, RK %d\n",Rstate.Pres, dir[i],RK);

		    printf("RK %d, Rstate.Pres <= 0, it's %g, 1 neighb dir %hd\n",RK, Rstate.Pres, dir[i]);
		    printf("Cell %hd (%e %e %e), type %hd, vol %g\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],
			   C->cell_type,C->cell_volume);
		    printf("Cell - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->Flow_data.State_vec.Rhof,
			   C->Flow_data.State_vec.rhof_prod,C->Flow_data.State_vec.RhoU,C->Flow_data.State_vec.RhoV,
			   C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE,C->Flow_data.State_vec.Pres);
		    printf("Cell prev - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->Flow_data.Temp_state_vec.Rhof,
			   C->Flow_data.Temp_state_vec.rhof_prod,C->Flow_data.Temp_state_vec.RhoU,C->Flow_data.Temp_state_vec.RhoV,
			   C->Flow_data.Temp_state_vec.RhoW,C->Flow_data.Temp_state_vec.RhoE,C->Flow_data.Temp_state_vec.Pres);
		    printf("Neigb - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Rhof,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.rhof_prod,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoU,C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoV,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoW,C->face_neighbours[dir[i]][0]->Flow_data.State_vec.RhoE,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Pres);
		    printf("Neigb prev - r %g, rp %g, ru %g, rv %g, rw %g, re %g, p %g\n",C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.Rhof,
			   C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.rhof_prod,
			   C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoU,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoV,
			   C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoW,C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.RhoE,
			   C->face_neighbours[dir[i]][0]->Flow_data.Temp_state_vec.Pres);
		    printf("Neighbour L %hd (%e %e %e), type %hd, vol %g, Pres %g\n",C->face_neighbours[dir[i]][0]->cell_level,
			   C->face_neighbours[dir[i]][0]->centroid[0],C->face_neighbours[dir[i]][0]->centroid[1],
			   C->face_neighbours[dir[i]][0]->centroid[2],
			   C->face_neighbours[dir[i]][0]->cell_type,C->face_neighbours[dir[i]][0]->cell_volume,
			   C->face_neighbours[dir[i]][0]->Flow_data.State_vec.Pres);
		    printf("%hd flx areas - %g %g %g %g\n",dir[i],C->flux_area[dir[i]][0],C->flux_area[dir[i]][1],C->flux_area[dir[i]][2],
			   C->flux_area[dir[i]][3]);
		    printf("opp_dir %hd flx areas - %g %g %g %g\n",opp_dir,C->face_neighbours[dir[i]][0]->flux_area[opp_dir][0],
			   C->face_neighbours[dir[i]][0]->flux_area[opp_dir][1],C->face_neighbours[dir[i]][0]->flux_area[opp_dir][2],
			   C->face_neighbours[dir[i]][0]->flux_area[opp_dir][3]);
		    
		    printf("Rstate - r %g, ru %g, rv %g, rw %g, re %g, rp %g, p %g\n",Rstate.Rhof,Rstate.RhoU,Rstate.RhoV,Rstate.RhoW,
			   Rstate.RhoE, Rstate.rhof_prod, Rstate.Pres);
		    printf("w_lim of neighb - %g, grad_w %g %g %g\n",C->face_neighbours[dir[i]][0]->Flow_data.Lim.w_lim,
			   C->face_neighbours[dir[i]][0]->Flow_data.Grads.grad_w[0],C->face_neighbours[dir[i]][0]->Flow_data.Grads.grad_w[1],
			   C->face_neighbours[dir[i]][0]->Flow_data.Grads.grad_w[2]);
		    if(C->Merge_data != NULL) {
		      printf("This cell's merged\n");
		      if(C->Merge_data == C->face_neighbours[dir[i]][0]->Merge_data)
			printf("Merged with neighb\n");
		    }

		    if(C->face_neighbours[dir[i]][0]->Merge_data != NULL)
		      printf("Neighbour's also merged\n");
		    if(C->is_small_cell == TRUE)
		      printf("C's a small cell\n");
		    if(C->face_neighbours[dir[i]][0]->is_small_cell == TRUE)
		      printf("Neighb's a small cell\n");

		    return(ERROR);
		  }
#endif

		  if(flux_type == ADAPTIVE) {
		    if((C -> shocked_cell == 'y') || (C->face_neighbours[dir[i]][0]->face_neighbours[dir[i]][0] == NULL))
		      compute_fluxes(&Lstate, &Rstate, norm, C -> Flow_data.Face_fluxes[dir[i]][0], EFM);
		    else compute_fluxes(&Lstate, &Rstate, norm, C -> Flow_data.Face_fluxes[dir[i]][0], AUSMDV);
		  }
		  else { /*Even if in that direction the neighbour is a border cell, use EFM*/
                    if(C->face_neighbours[dir[i]][0]->face_neighbours[dir[i]][0] == NULL)
                      compute_fluxes(&Lstate, &Rstate, norm, C -> Flow_data.Face_fluxes[dir[i]][0], EFM);
                    else compute_fluxes(&Lstate, &Rstate, norm, C -> Flow_data.Face_fluxes[dir[i]][0], flux_type);
		  }
		}	  
	      else if(C -> face_neighbours[dir[i]][1] != NULL) /*4 cells exist in this direction*/
		{
		  reconstructed_yet = FALSE;

		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Cycle through each quadrant*/
		    {
		      if(
#if DETON
			 (C -> face_neighbours[dir[i]][j] -> un_det != TRUE) &&
#endif
			 (C -> Flow_data.Face_fluxes[dir[i]][j] != NULL)) /*Flux vector already exists*/
			{
			  if(reconstructed_yet == FALSE)
			    {
			      Lstate = reconstruct_flow(C, dir[i], j, IRRELEVANT, IRRELEVANT, EOF, NULL, NULL, recon_full);
			      reconstructed_yet = TRUE;
			    }

			  Rstate = reconstruct_flow(C -> face_neighbours[dir[i]][j], opp_dir, IRRELEVANT, IRRELEVANT, IRRELEVANT, EOF, 
						    NULL, NULL, 'n');

#if 0
			  if(((Rstate.RhoE) - 0.5*(SQR(Rstate.RhoU) + SQR(Rstate.RhoV) + SQR(Rstate.RhoW))/(Rstate.Rhof)) < 0) {
			    printf("It's here 4, C %hd (%9.8e %9.8e %9.8e), dir %hd, neighb no. %hd (%9.8e %9.8e %9.8e)\n",C->cell_level,C->centroid[0],
				   C->centroid[1],C->centroid[2],dir[i],j,C->face_neighbours[dir[i]][j]->centroid[0],
				   C->face_neighbours[dir[i]][j]->centroid[1],C->face_neighbours[dir[i]][j]->centroid[2]);
			    printf("C's vol %g, type %hd, un_det %hd, neighb's vol %g, type %hd, un_det %hd\n",C->cell_volume,C->cell_type,C->un_det,C->face_neighbours[dir[i]][j]->cell_volume,C->face_neighbours[dir[i]][j]->cell_type,C->face_neighbours[dir[i]][j]->un_det);
			    printf("Time is %e, stage is %hd, explosion radisu should be %g\n",time,stage,6717.4*time);
			    printf("C's r %g, ru %g, rv %g, rw %g, re %g, rp %g\n",C->Flow_data.State_vec.Rhof,C->Flow_data.State_vec.RhoU,
				   C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoW,C->Flow_data.State_vec.RhoE,
				   C->Flow_data.State_vec.rhof_prod);
			    printf("Lstate r %g, ru %g, rv %g, rw %g, re %g, rp %g\n",Lstate.Rhof,Lstate.RhoU,Lstate.RhoV,Lstate.RhoW,
				   Lstate.RhoE,Lstate.rhof_prod);
			    printf("Rstate -> r %g, ru %g, rv %g, rw %g, re %g, rp %g\n",Rstate.Rhof,Rstate.RhoU,Rstate.RhoV,Rstate.RhoW,
				   Rstate.RhoE,Rstate.rhof_prod);
			    printf("Neighb's r %g, ru %g, rv %g, rw %g, re %g, rp %g\n",C->face_neighbours[dir[i]][j]->Flow_data.State_vec.Rhof,C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoU,C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoV,C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoW,C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoE,C->face_neighbours[dir[i]][j]->Flow_data.State_vec.rhof_prod);
			    return(ERROR);
			  }
#endif

#if 0
			  if(Lstate.Rhof <= 0) { 
			    printf("Lstate.Rhof <= 0, it's %g, 4 neighb dir %hd, quad %hd\n",Lstate.Rhof, dir[i], j);
			    return(ERROR);
			  }
			  else if(Lstate.Rhof > 0) {
			  }
			  else{
			    printf("Lstate Rhof bad - %g, dir %hd, quad %hd, RK %d\n",Lstate.Rhof,dir[i],j,RK);
			    return(ERROR);
			  }

			  if(Rstate.Rhof <= 0) {
			    printf("Rstate.Rhof <= 0, it's %g, 4 neighb dir %hd, quad %hd, type %hd\n",Rstate.Rhof, dir[i], j,
				   C->face_neighbours[dir[i]][j]->cell_type);
			    printf("C %hd, (%e %e %e), type %hd, flx_area[%hd][%hd] %g\n",C->cell_level,C->centroid[0],C->centroid[1],
				   C->centroid[2],C->cell_type,dir[i],j,C->flux_area[dir[i]][j]);			
			    printf("Neighb %hd (%e %e %e),flx_area %g\n",C->face_neighbours[dir[i]][j]->cell_level,
				   C->face_neighbours[dir[i]][j]->centroid[0],C->face_neighbours[dir[i]][j]->centroid[1],
				   C->face_neighbours[dir[i]][j]->centroid[2],
				   (C->face_neighbours[dir[i]][j]->flux_area[opp_dir][0]+
				    C->face_neighbours[dir[i]][j]->flux_area[opp_dir][1]+C->face_neighbours[dir[i]][j]->flux_area[opp_dir][2]+
				    C->face_neighbours[dir[i]][j]->flux_area[opp_dir][2]));
			    printf("Neighb par %hd, (%e %e %e)\n",C->face_neighbours[dir[i]][j]->parent->cell_level,
				   C->face_neighbours[dir[i]][j]->parent->centroid[0],C->face_neighbours[dir[i]][j]->parent->centroid[1],
				   C->face_neighbours[dir[i]][j]->parent->centroid[2]);
			    printf("Neighb r %g, p %g, ru %g, rv %g, rw %g, re %g\n",C->face_neighbours[dir[i]][j]->Flow_data.State_vec.Rhof,
				   C->face_neighbours[dir[i]][j]->Flow_data.State_vec.Pres,
				   C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoU,
				   C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoV,
				   C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoW,
				   C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoE);
			    return(ERROR);
			  }
			  else if(Rstate.Rhof > 0) {
			  }
			  else {
			    printf("Rstate Rhof bad - %g, dir %hd, quad %hd, RK %d\n",Rstate.Rhof, dir[i],j,RK);			
			    return(ERROR);
			  }

			  if(Rstate.Pres <= 0) {
			    printf("RK %d\n", RK);
			    printf("Rstate.Pres <= 0, it's %g, 4 neighb dir %hd, quad %hd\n",Rstate.Pres, dir[i], j);
			    printf("Rstate.Rhof %g\n",Rstate.Rhof);
			    printf("Rstate.RhoE %g\n",Rstate.RhoE);
			    printf("Rstate.RhoU %g\n",Rstate.RhoU);
			    printf("Rstate.RhoV %g\n",Rstate.RhoV);
			    printf("Rstate.RhoW %g\n",Rstate.RhoW);

			    printf("Cell type is level %hd, (%g %g %g), type %hd, volume %g, small_cell status %hd\n",C->cell_level,C->centroid[0],
				   C->centroid[1],C->centroid[2],C->cell_type,C->cell_volume,C->is_small_cell);

			    printf("Cell's flux area dir %hd quad %hd is %g\n",i,j,C->flux_area[dir[i]][j]);

			    printf("Cell's p is %g\n",C->Flow_data.State_vec.Pres);
			    printf("Cell's r is %g\n",C->Flow_data.State_vec.Rhof);
			    printf("Cell's ru is %g\n",C->Flow_data.State_vec.RhoU);
			    printf("Cell's rv is %g\n",C->Flow_data.State_vec.RhoV);
			    printf("Cell's rw is %g\n",C->Flow_data.State_vec.RhoW);
			    printf("Cell's re is %g\n",C->Flow_data.State_vec.RhoE);

			    printf("Cell's parent level %hd, (%g %g %g), type %hd\n",C->parent->cell_level,
				   C->parent->centroid[0],C->parent->centroid[1],C->parent->centroid[2],C->parent->cell_type);

			    if(C->Merge_data != NULL)
			      printf("This cell is merged\n");

			    if(C->face_neighbours[dir[i]][j] == NULL)
			      printf("Cell's neighbour dir %hd quad %hd is NULL\n",dir[i],j);
			    else printf("Neighbour in dir %hd quad %hd does exist\n",dir[i],j);
		
			    if(C->face_neighbours[dir[i]][0] == NULL)
			      printf("Cell is on border in this dir %hd\n",dir[i]);
		
			    if(C->face_neighbours[dir[i]][1] == NULL)
			      printf("No neighb dir %hd quad 1\n",dir[i]);

			    if(C->face_neighbours[dir[i]][2] == NULL)
			      printf("No neighb dir %hd quad 2\n",dir[i]);

			    printf("Neighbour's level %hd, (%g %g %g), type %hd, volume %g, small_cell_status %hd\n",
				   C->face_neighbours[dir[i]][j]->cell_level,C->face_neighbours[dir[i]][j]->centroid[0],C->face_neighbours[dir[i]][j]->centroid[1],
				   C->face_neighbours[dir[i]][j]->centroid[2],C->face_neighbours[dir[i]][j]->cell_type,C->face_neighbours[dir[i]][j]->cell_volume,
				   C->face_neighbours[dir[i]][j]->is_small_cell);

			    printf("Neighbour's parent level %hd, (%g %g %g), type %hd\n",C->face_neighbours[dir[i]][j]->parent->cell_level,
				   C->face_neighbours[dir[i]][j]->parent->centroid[0],C->face_neighbours[dir[i]][j]->parent->centroid[1],
				   C->face_neighbours[dir[i]][j]->parent->centroid[2],C->face_neighbours[dir[i]][j]->parent->cell_type);

			    printf("Neighbour's p is %g\n",C->face_neighbours[dir[i]][j]->Flow_data.State_vec.Pres);
			    printf("Neighbour's r is %g\n",C->face_neighbours[dir[i]][j]->Flow_data.State_vec.Rhof);
			    printf("Neighbour's ru is %g\n",C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoU);
			    printf("Neighbour's rv is %g\n",C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoV);
			    printf("Neighbour's rw is %g\n",C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoW);
			    printf("Neighbour's re is %g\n",C->face_neighbours[dir[i]][j]->Flow_data.State_vec.RhoE);

			    printf("Neighbour's flux area in opposite dir %hd is %g\n",opp_dir,
				   ((C->face_neighbours[dir[i]][j]->flux_area[opp_dir][0]) +
				    (C->face_neighbours[dir[i]][j]->flux_area[opp_dir][1]) + (C->face_neighbours[dir[i]][j]->flux_area[opp_dir][2]) +
				    (C->face_neighbours[dir[i]][j]->flux_area[opp_dir][3])));

			    if(C->face_neighbours[dir[i]][j] != NULL)
			      printf("This neighbour is merged\n");

			    if(C->Merge_data == C->face_neighbours[dir[i]][j]->Merge_data) 
			      printf("They're merged to each other");
			    else printf("They aren't merged to each other\n");		

			    return(ERROR);
			  }
			  else if(Rstate.Pres > 0) {
			  }
			  else {
			    printf("Rstate Pres bad - %g, dir %hd, quad %hd, RK %d\n",Rstate.Pres, dir[i],j,RK);
			    return(ERROR);
			  }
#endif

			  if(flux_type == ADAPTIVE) {
			    if((C -> shocked_cell == 'y') || (C->face_neighbours[dir[i]][j]->face_neighbours[dir[i]][0] == NULL))
			      compute_fluxes(&Lstate, &Rstate, norm, C -> Flow_data.Face_fluxes[dir[i]][j], EFM);
			    else compute_fluxes(&Lstate, &Rstate, norm, C -> Flow_data.Face_fluxes[dir[i]][j], AUSMDV);
			  }
			  else { /*Even if in that direction the neighbour is a border cell, use EFM*/
			    if(C->face_neighbours[dir[i]][j]->face_neighbours[dir[i]][0] == NULL)
			      compute_fluxes(&Lstate, &Rstate, norm, C -> Flow_data.Face_fluxes[dir[i]][j], EFM);
			    else compute_fluxes(&Lstate, &Rstate, norm, C -> Flow_data.Face_fluxes[dir[i]][j], flux_type);
			  }
			} 		  
		    }
		} 
	    }
#if DETON 
	}
#endif

      if(node == Tail)
	break;
      else node = node -> next;      
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Compute fluxes between 2 cells with their corresp. state vectors at the (generalized) interface.  norm is outward normal
   from left face*/

void compute_fluxes(State_vector * Lstate, State_vector * Rstate, double norm[], Flux_vector flux_values, short int type)
{   
  if(type == AUSMDV) {
    double Pl, Pr, Ul[3], Ur[3], Vl[3], Vr[3], ur, ul, ulp, urm, Al, Ar, rl, rr, al, ar, am, d1, Hl, Hr, Hdash, Gdash, G_pdash;
    double Plp, Prm, Phalf, Lddash, Lvdash, Lndash, Ltdash[3], Ldash[3], s, K, C, delua;

#if AFTERBURN
    double G_pubdash;
#endif
    G_pdash = 0;
    K = 10.0;

    rl = Lstate->Rhof;
    rr = Rstate->Rhof;

    Ul[0] = (Lstate->RhoU)/(rl);
    Ul[1] = (Lstate->RhoV)/(rl);
    Ul[2] = (Lstate->RhoW)/(rl);

    Ur[0] = (Rstate->RhoU)/(rr);
    Ur[1] = (Rstate->RhoV)/(rr);
    Ur[2] = (Rstate->RhoW)/(rr);

    Pl = Lstate -> Pres; 
    Pr = Rstate -> Pres;

    if(gas_data.num_species > 1)
      {      
	al = get_mixture_soundspd(Lstate);
	ar = get_mixture_soundspd(Rstate);
      }
    else
      { 
	al = SQRT((gas_data.gam_amb)*Pl/(rl)); 
	ar = SQRT((gas_data.gam_amb)*Pr/(rr));
      }

#if 0
    if(al <= 0)
      printf("We gotta bad soundspd l, %g, Pl %g, Rhofl %g\n",al,Pl,Lstate->Rhof);
    else if(al > 0) {
    }
    else {
      printf("We gotta bad soundspdl, %g, Pl %g, Rhofl %g\n",al,Pl,Lstate->Rhof);
      abort();
    }

    if(ar <= 0)
      printf("We gotta bad soundspd r, %g, Pr %g, Rhofr %g\n",ar,Pr,Rstate->Rhof);
    else if(ar > 0) {
    }
    else {
      printf("We gotta bad soundspdr, %g, Pr %g, Rhofr %g\n",ar,Pr,Rstate->Rhof);
      abort();
    }
#endif

    ul = Ul[0]*norm[0] + Ul[1]*norm[1] + Ul[2]*norm[2];
    ur = Ur[0]*norm[0] + Ur[1]*norm[1] + Ur[2]*norm[2];

    Vl[0] = Ul[0] - ul*norm[0];
    Vl[1] = Ul[1] - ul*norm[1];
    Vl[2] = Ul[2] - ul*norm[2];

    Vr[0] = Ur[0] - ur*norm[0];
    Vr[1] = Ur[1] - ur*norm[1];
    Vr[2] = Ur[2] - ur*norm[2];

    d1 = (Pl/rl) + (Pr/rr);
    Al = (2.0*Pl/rl)/d1; 
    Ar = (2.0*Pr/rr)/d1;

    am = MAX(al,ar);

    d1 = (ul + ABS(ul))/2.0;
    if((ABS(ul)/am) <= 1.0) 
      ulp = Al*(SQR(ul+am)/(4.0*am) - d1) + d1;
    else ulp = d1;
  
    d1 = (ur - ABS(ur))/2.0;
    if((ABS(ur)/am) <= 1.0) 
      urm = -Ar*(SQR(ur-am)/(4.0*am) + d1) + d1;
    else urm = d1;

    Hl = ((Lstate->RhoE) + Pl)/rl;
    Hr = ((Rstate->RhoE) + Pr)/rr;

    Gdash = ulp*rl + urm*rr; /*Mass flux*/
    if(gas_data.num_species > 1) {
      G_pdash = ulp*(Lstate->rhof_prod) + urm*(Rstate->rhof_prod);

#if AFTERBURN
      G_pubdash = ulp*(Lstate->rhof_prod_ub) + urm*(Rstate->rhof_prod_ub);
#endif
    }
  
    Hdash = 0.5*(Gdash*(Hl+Hr) - ABS(Gdash)*(Hr-Hl)); /*Enthalpy flux*/

#if 1
    /*2nd order pressure splitting*/
  
    if((ABS(ul)/am) <= 1.0) 
      Plp = 0.25*Pl*SQR((ul/am)+1.0)*(2.0-(ul/am));
    else Plp = Pl*(ul + ABS(ul))/(2.0*ul);

    if((ABS(ur)/am) <= 1.0) 
      Prm = 0.25*Pr*SQR((ur/am)-1.0)*(2.0+(ur/am));
    else Prm = Pr*(ur - ABS(ur))/(2.0*ur);
  
#else
    /*1st order pressure splitting*/

    if((ABS(ul)/am) <= 1.0)
      Plp = 0.5*Pl*(1.0+ul/am);
    else Plp = Pl*(ul+ABS(ul))/(2.0*ul);

    if((ABS(ur)/am) <= 1.0)
      Prm = 0.5*Pr*(1.0-ur/am);
    else Prm = Pr*(ur-ABS(ur))/(2.0*ur);
#endif

    Phalf = Plp + Prm;

    s = 0.5*MIN(1.0, K*ABS(Pr-Pl)/MIN(Pl,Pr)); /*Switching factor*/

    Lddash = 0.5*(Gdash*(ul+ur) - ABS(Gdash)*(ur-ul));
    Lvdash = rl*ul*ulp + rr*ur*urm;
    Lndash = (0.5+s)*Lvdash + (0.5-s)*Lddash + Phalf;
  
    Ltdash[0] = 0.5*(Gdash*(Vl[0]+Vr[0]) - ABS(Gdash)*(Vr[0]-Vl[0]));
    Ltdash[1] = 0.5*(Gdash*(Vl[1]+Vr[1]) - ABS(Gdash)*(Vr[1]-Vl[1]));
    Ltdash[2] = 0.5*(Gdash*(Vl[2]+Vr[2]) - ABS(Gdash)*(Vr[2]-Vl[2]));
  
    Ldash[0] = Lndash*norm[0] + Ltdash[0];
    Ldash[1] = Lndash*norm[1] + Ltdash[1];
    Ldash[2] = Lndash*norm[2] + Ltdash[2]; /*Momentum flux*/

    flux_values -> Mss_flx = Gdash;
    flux_values -> X_mntm_flx = Ldash[0];
    flux_values -> Y_mntm_flx = Ldash[1];
    flux_values -> Z_mntm_flx = Ldash[2];
    flux_values -> E_flx = Hdash;

    if(gas_data.num_species > 1) {
      flux_values -> mss_flx_prod = G_pdash;

#if AFTERBURN
      flux_values -> mss_flx_prod_ub = G_pubdash;
#endif
    }

    /*Now to fix the glitch at the expansion sonic point*/

    if(((ul - al) < 0) && ((ur - ar) > 0)) {
      C = 0.125;
      delua = (ur - ar) - (ul - al);
      flux_values -> Mss_flx = Gdash - C*delua*(rr - rl);
      flux_values -> X_mntm_flx = Ldash[0] - C*delua*((Rstate->RhoU) - (Lstate->RhoU));
      flux_values -> Y_mntm_flx = Ldash[1] - C*delua*((Rstate->RhoV) - (Lstate->RhoV));
      flux_values -> Z_mntm_flx = Ldash[2] - C*delua*((Rstate->RhoW) - (Lstate->RhoW));
      flux_values -> E_flx = Hdash - C*delua*(((Rstate->RhoE) + Pr) - ((Lstate->RhoE) + Pl));
      if(gas_data.num_species > 1) {
	flux_values -> mss_flx_prod = G_pdash - C*delua*((Rstate->rhof_prod) - (Lstate->rhof_prod));
#if AFTERBURN
	flux_values -> mss_flx_prod_ub = G_pubdash - C*delua*((Rstate->rhof_prod_ub) - (Rstate->rhof_prod_ub));
#endif
      } 
    }
    else if(((ul + al) < 0) && ((ur + ar) > 0)) {
      C = 0.125;
      delua = (ur + ar) - (ul + al);
      flux_values -> Mss_flx = Gdash - C*delua*(rr - rl);
      flux_values -> X_mntm_flx = Ldash[0] - C*delua*((Rstate->RhoU) - (Lstate->RhoU));
      flux_values -> Y_mntm_flx = Ldash[1] - C*delua*((Rstate->RhoV) - (Lstate->RhoV));
      flux_values -> Z_mntm_flx = Ldash[2] - C*delua*((Rstate->RhoW) - (Lstate->RhoW));
      flux_values -> E_flx = Hdash - C*delua*(((Rstate->RhoE) + Pr) - ((Lstate->RhoE) + Pl));
      if(gas_data.num_species > 1) {
	flux_values -> mss_flx_prod = G_pdash - C*delua*((Rstate->rhof_prod) - (Lstate->rhof_prod));
#if AFTERBURN
	flux_values -> mss_flx_prod_ub = G_pubdash - C*delua*((Rstate->rhof_prod_ub) - (Rstate->rhof_prod_ub));
#endif
      }
    }
  }
  else if(type == AUSM) {
    double Pl, Pr, soundspdl, soundspdr, Ml, Mr, Mhalf, Mlplus, Mrminus, Plplus, Prminus, left_coeff, right_coeff, Ul[3], Ur[3];

    Ul[0] = (Lstate->RhoU)/(Lstate->Rhof);
    Ul[1] = (Lstate->RhoV)/(Lstate->Rhof);
    Ul[2] = (Lstate->RhoW)/(Lstate->Rhof);

    Ur[0] = (Rstate->RhoU)/(Rstate->Rhof);
    Ur[1] = (Rstate->RhoV)/(Rstate->Rhof);
    Ur[2] = (Rstate->RhoW)/(Rstate->Rhof);

    Pl = Lstate -> Pres; 
    Pr = Rstate -> Pres;

    if(gas_data.num_species > 1)
      {      
	soundspdl = get_mixture_soundspd(Lstate);
	soundspdr = get_mixture_soundspd(Rstate);
      }
    else
      { 
	soundspdl = SQRT((gas_data.gam_amb)*Pl/(Lstate -> Rhof)); 
	soundspdr = SQRT((gas_data.gam_amb)*Pr/(Rstate -> Rhof));
      }

    Ml = (Ul[0]*norm[0] + Ul[1]*norm[1] + Ul[2]*norm[2])/soundspdl;
    Mr = (Ur[0]*norm[0] + Ur[1]*norm[1] + Ur[2]*norm[2])/soundspdr;

#if 0
    if(soundspdl <= 0)
      printf("We gotta bad soundspd l, %g, Pl %g, Rhofl %g\n",soundspdl,Pl,Lstate->Rhof);
    else if(soundspdl > 0) {
    }
    else {
      printf("We gotta bad soundspdl, %g, Pl %g, Rhofl %g\n",soundspdl,Pl,Lstate->Rhof);
    }

    if(soundspdr <= 0)
      printf("We gotta bad soundspd r, %g, Pr %g, Rhofr %g\n",soundspdr,Pr,Rstate->Rhof);
    else if(soundspdr > 0) {
    }
    else {
      printf("We gotta bad soundspdr, %g, Pr %g, Rhofr %g\n",soundspdr,Pr,Rstate->Rhof);
    }
#endif
  
    if(ABS(Ml) <= 1.0)
      {
	Mlplus = 0.25*SQR(Ml + 1.0);
	Plplus = 0.5*Pl*(1.0 + Ml);
      }
    else 
      {
	Mlplus = 0.5*(Ml + ABS(Ml)); 
	Plplus = 0.5*Pl*(1.0 + ABS(Ml)/Ml);
      }

    if(ABS(Mr) <= 1.0)
      {
	Mrminus = -0.25*SQR(Mr - 1.0);
	Prminus = 0.5*Pr*(1.0 - Mr);
      }
    else
      { 
	Mrminus = 0.5*(Mr - ABS(Mr)); 
	Prminus = 0.5*Pr*(1.0 - ABS(Mr)/Mr);
      }
  
    Mhalf = Mlplus + Mrminus;
    left_coeff = 0.5*(Mhalf + ABS(Mhalf))*soundspdl;
    right_coeff = 0.5*(Mhalf - ABS(Mhalf))*soundspdr;

    flux_values -> Mss_flx = left_coeff*(Lstate->Rhof) + right_coeff*(Rstate->Rhof);
    flux_values -> X_mntm_flx = left_coeff*(Lstate->RhoU) + right_coeff*(Rstate->RhoU) + norm[0]*(Plplus + Prminus);
    flux_values -> Y_mntm_flx = left_coeff*(Lstate->RhoV) + right_coeff*(Rstate->RhoV) + norm[1]*(Plplus + Prminus);
    flux_values -> Z_mntm_flx = left_coeff*(Lstate->RhoW) + right_coeff*(Rstate->RhoW) + norm[2]*(Plplus + Prminus);
    flux_values -> E_flx = left_coeff*((Lstate->RhoE) + Pl) + right_coeff*((Rstate->RhoE) + Pr);

    if(gas_data.num_species > 1) {
      flux_values -> mss_flx_prod = left_coeff*(Lstate->rhof_prod) + right_coeff*(Rstate->rhof_prod); /*Species (prodcuts gas) continuity*/
#if AFTERBURN
      flux_values -> mss_flx_prod_ub = left_coeff*(Lstate->rhof_prod_ub) + right_coeff*(Rstate->rhof_prod_ub);
#endif
    } 

#if 0
    if(flux_values -> Mss_flx < 0) {
#if 0
      printf("Saw flux_values -> Mss_flx %g\n",flux_values->Mss_flx);

      if((direction == NTH) || (direction == STH)) /*Y direction*/
	{
	  printf("R.RhoV is %g, R.Rhof is %g\n",Rstate->RhoV,Rstate->Rhof);
	}
      else if((direction == EST) || (direction == WST)) /*X direction*/
	{
	  printf("R.RhoU is %g, R.Rhof is %g\n",Rstate->RhoU,Rstate->Rhof);
	}
      else /*Z direction*/
	{
	  printf("R.RhoW is %g, R.Rhof is %g\n",Rstate->RhoW,Rstate->Rhof);
	  if(Rstate->RhoW <= 0) {
	    printf("It's <= 0 ");
	    if(Rstate->RhoW <= -DBL_EPSILON)
	      printf("and also too small\n");
	    else printf("\n");
	  }
	  else if(Rstate->RhoW > 0)
	    printf("It's > 0\n");
	  else printf("It's truly nan\n");
	}
#endif	 
    }
    else if(flux_values -> Mss_flx >= 0) {
    }
    else {
#if 0
      printf("Saw flux_values -> Mss_flx %g\n",flux_values->Mss_flx);
#endif
    }

    if(flux_values -> E_flx < 0)
      {}
    else if(flux_values -> E_flx >= 0) {
    }
    else {
      printf("Saw flux_values -> E_flx %g\n",flux_values->E_flx);
    }

    if(flux_values -> X_mntm_flx < 0)
      {}
    else if(flux_values -> X_mntm_flx >= 0) {
    }
    else {
      printf("Saw flux_values -> X_mntm_flx %g\n",flux_values->X_mntm_flx);
    }
 
    if(flux_values -> Y_mntm_flx < 0)
      {}
    else if(flux_values -> Y_mntm_flx >= 0) {
    }
    else {
      printf("Saw flux_values -> Y_mntm_flx %g\n",flux_values->Y_mntm_flx);
    }

    if(flux_values -> Z_mntm_flx < 0)
      {}
    else if(flux_values -> Z_mntm_flx >= 0) {
    }
    else {
      printf("Saw flux_values -> Z_mntm_flx %g\n",flux_values->Z_mntm_flx);
    }
#endif
  }
  else if(type == EFM) {
    /*Pretty much a copy of Paul Petrie's code except it's adapted for generalized interfaces*/
  
    double uL, vL, uR, vR, wL, wR, VLs, VRs, vnL, vnR, tL, tR;
    double rtL, cmpL, rtR, cmpR, R_amb;
    double hvsqL, hvsqR;
    double WL, WR, dL, dR;
    double rhoL, rhoR, presL, presR;
    double hL, hR;
    double snL, snR, exL, exR, efL, efR;
    double fmsL, fmsR;
    double cv, cp, con, gam, Rgas;
    double cvL, cvR, RgasL, RgasR;
    double rLsqrt, rRsqrt, alpha;
    double pmfL, pmfR;

    /*Calculate Constants*/
    /* dtwspi = 1.0 / (2.0 * sqrt ( PI )); */
#   define  dtwspi  0.282094792
#   define PHI 1.0
    pmfL = 0; pmfR = 0;
    rhoL = Lstate -> Rhof;
    presL = Lstate -> Pres;
    uL = (Lstate -> RhoU)/rhoL;
    vL = (Lstate -> RhoV)/rhoL;
    wL = (Lstate -> RhoW)/rhoL;
    VLs = SQR(uL) + SQR(vL) + SQR(wL);
    hL = ((Lstate -> RhoE)/rhoL) - 0.5*VLs + (presL/rhoL);

    rhoR = Rstate -> Rhof;
    presR = Rstate -> Pres;
    uR = (Rstate -> RhoU)/rhoR;
    vR = (Rstate -> RhoV)/rhoR;
    wR = (Rstate -> RhoW)/rhoR;
    VRs = SQR(uR) + SQR(vR) + SQR(wR);
    hR = ((Rstate -> RhoE)/rhoR) - 0.5*VRs + (presR/rhoR);

    /* Derive the gas "constants" from the local conditions. */
    if(gas_data.num_species > 1) {
      pmfL = (Lstate -> rhof_prod)/rhoL;
      pmfR = (Rstate -> rhof_prod)/rhoR;
      cvL = (1.0-pmfL)*(gas_data.Cv_amb) + pmfL*(gas_data.Cv_prod);
      cvR = (1.0-pmfL)*(gas_data.Cv_amb) + pmfL*(gas_data.Cv_prod);
      tL = get_mixture_T(rhoL, Lstate -> rhof_prod, Lstate -> RhoE, presL, uL, vL, wL);
      tR = get_mixture_T(rhoR, Rstate -> rhof_prod, Rstate -> RhoE, presR, uR, vR, wR);
    }
    else { /*Same 'ambient' species*/
      cvL = gas_data.Cv_amb;
      cvR = cvL;
      R_amb = (gas_data.Cv_amb)*((gas_data.gam_amb) - 1.0); /*Ambient gas R*/
      tL = presL/(R_amb*rhoL);
      tR = presR/(R_amb*rhoR);
    }

    RgasL = presL / (rhoL * tL); /*Gas 'constant' for left state gas*/
    RgasR = presR / (rhoR * tR);
                                                                                                                          
    rLsqrt = sqrt(rhoL);
    rRsqrt = sqrt(rhoR);
    alpha = rLsqrt / (rLsqrt + rRsqrt);
                                                                                                                          
    cv = alpha * cvL + (1.0 - alpha) * cvR;
    Rgas = alpha * RgasL + (1.0 - alpha) * RgasR;
                                                                                                                          
    cp = cv + Rgas;
    gam = cp / cv;

    /*Start EFM calculation proper*/
    vnL = uL*norm[0] + vL*norm[1] + wL*norm[2]; /*Get velocity normal to interface*/
    vnR = uR*norm[0] + vR*norm[1] + wR*norm[2];

    con = 0.5 * (gam + 1.0) / (gam - 1.0);

#if 0  
    rtL = Rgas * tL; /*This is how Paul Petrie coded it, but shouldn't it be presL/rhoL (like with rtR)?*/
#else
    rtL = presL / rhoL;
#endif
    cmpL = sqrt(2.0 * rtL);
    hvsqL = 0.5 * VLs;
  
    snL = vnL / (PHI * cmpL);
    exxef(snL, &exL, &efL);
  
    WL = 0.5 * (1.0 + efL);
    dL = exL * dtwspi;
  
    rtR = presR / rhoR;
    cmpR = sqrt(2.0 * rtR);
    hvsqR = 0.5 * VRs;
  
    snR = vnR / (PHI * cmpR);
    exxef(snR, &exR, &efR);
  
    WR = 0.5 * (1.0 - efR);
    dR = -exR * dtwspi;

    /*Combine the fluxes*/
    fmsL = (WL * rhoL * vnL) + (dL * cmpL * rhoL);
    fmsR = (WR * rhoR * vnR) + (dR * cmpR * rhoR);
  
    flux_values -> Mss_flx = fmsL + fmsR;
  
    flux_values -> X_mntm_flx = fmsL * uL + fmsR * uR + norm[0]*(WL * presL + WR * presR);
    flux_values -> Y_mntm_flx = fmsL * vL + fmsR * vR + norm[1]*(WL * presL + WR * presR);
    flux_values -> Z_mntm_flx = fmsL * wL + fmsR * wR + norm[2]*(WL * presL + WR * presR);
  
    flux_values -> E_flx = (WL * rhoL * vnL) * (hvsqL + hL) +
      (WR * rhoR * vnR) * (hvsqR + hR) +
      (dL * cmpL * rhoL) * (hvsqL + con * rtL) +
      (dR * cmpR * rhoR) * (hvsqR + con * rtR);
  
    if(gas_data.num_species > 1) {
      if ((flux_values -> Mss_flx) > 0.0)
	flux_values -> mss_flx_prod = (flux_values -> Mss_flx) * pmfL;
      else flux_values -> mss_flx_prod = (flux_values -> Mss_flx) * pmfR;
    }
  }
}

/*------------------------------------------------------------------*/

/**\brief Compute minimum local timestep for a cell by comparing time taken for fastest wave to cross
   a cell (assuming CFL of 1) - only 'normal' cells will execute this function*/

double min_local_cell_timestep(Cart_cell C)
{
  short int i, j;
  double Vel, soundspd, sum_dem, Vol, dx, shock_param;
  
  /*We reconstruct velocities to cell interfaces (don't worry if the interface shared by > 2 cells as we're
    only looking at time to cross a single cell)*/

  sum_dem = 0; /*Sum over all interface areas multiplied by maximum wave speeds*/
  shock_param = 0;
  if(gas_data.num_species > 1)
    soundspd = get_mixture_soundspd(&(C -> Flow_data.State_vec));
  else soundspd = SQRT((gas_data.gam_amb)*(C->Flow_data.State_vec.Pres)/(C->Flow_data.State_vec.Rhof));

  if(flux_type == ADAPTIVE) {

    /*Calculated 'shocked cells' to turn off reconstruction*/

    dx = CALC_CELL_EDGE_LENGTH(C);
    C -> shocked_cell = 'n';

    for(i = 0; i < 3; i++) {
      switch(i) {
      case 0: /*X direction*/
	shock_param = (C->Flow_data.Grads.grad_u[0])*dx/soundspd; 
	break;
      case 1:
	shock_param = (C->Flow_data.Grads.grad_v[1])*dx/soundspd;
	break;
      case 2:
	shock_param = (C->Flow_data.Grads.grad_w[2])*dx/soundspd;
	break;
      }

      if(shock_param <= -adapt_param) {
	C -> shocked_cell = 'y';      
	break;
      }
    }
  }

  if(C -> Merge_data == NULL)
    {
      for(i = 0; i < NUM_FACES; i++)
	{
	  if((i == EST) || (i == WST))
	    reconstruct_flow(C, i, IRRELEVANT, IRRELEVANT, IRRELEVANT, 'U', &Vel, NULL, 'n'); /*Only want velocities*/
	  else if((i == NTH) || (i == STH))
	    reconstruct_flow(C, i, IRRELEVANT, IRRELEVANT, IRRELEVANT, 'V', &Vel, NULL, 'n');
	  else reconstruct_flow(C, i, IRRELEVANT, IRRELEVANT, IRRELEVANT, 'W', &Vel, NULL, 'n');
	 
	  sum_dem += (C->flux_area[i][0] + C->flux_area[i][1] + C->flux_area[i][2] + C->flux_area[i][3])*(soundspd + ABS(Vel));
	}

      Vol = C -> cell_volume;
    }
  else
    {
      List_merge L = C -> Merge_data -> Linked_cells;
          
      /*Go through list of linked cells and sum area-wave products over all non-shared interfaces*/
  
      while(L != NULL)
	{
	  for(i = 0; i < NUM_FACES; i++)
	    {
	      if((i == EST) || (i == WST))
		reconstruct_flow(L -> cell_loc, i, IRRELEVANT, IRRELEVANT, IRRELEVANT, 'U', &Vel, NULL, 'n');
	      else if((i == NTH) || (i == STH))
		reconstruct_flow(L -> cell_loc, i, IRRELEVANT, IRRELEVANT, IRRELEVANT, 'V', &Vel, NULL, 'n');
	      else reconstruct_flow(L -> cell_loc, i, IRRELEVANT, IRRELEVANT, IRRELEVANT, 'W', &Vel, NULL, 'n');

	      if(L -> cell_loc -> face_neighbours[i][0] == NULL)
		sum_dem += (L->cell_loc->flux_area[i][0] + L->cell_loc->flux_area[i][1] + L->cell_loc->flux_area[i][2] +
			    L->cell_loc->flux_area[i][3])*(soundspd + ABS(Vel));	
	      else if(L -> cell_loc -> face_neighbours[i][1] == NULL)
		{ 
		  if(L -> cell_loc -> face_neighbours[i][0] -> Merge_data != L -> cell_loc -> Merge_data) 
		    { /*This interface isn't shared by another linked cell*/

		      sum_dem += (L->cell_loc->flux_area[i][0] + L->cell_loc->flux_area[i][1] + L->cell_loc->flux_area[i][2] +
				  L->cell_loc->flux_area[i][3])*(soundspd + ABS(Vel));
		    } /*If the interface is shared, we shouldn't sum it*/
		}
	      else /*4 cells in this direction*/
		{
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		    {
		      if(L -> cell_loc -> face_neighbours[i][j] -> Merge_data != L -> cell_loc -> Merge_data)
			sum_dem += (L -> cell_loc -> flux_area[i][j])*(soundspd + ABS(Vel));
		    }
		}	      
	    }

	  L = L -> next;
	}

      Vol = C -> Merge_data -> cell_volume;
    }

#if 0
  if(Vol/sum_dem <= 0) {
    
    printf("Got %g time, C %hd (%12.11e %12.11e %12.11e), type %hd\n",Vol/sum_dem,C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],
	   C->cell_type);
    printf("Vol is %g, sum_dem is %g\n",Vol,sum_dem);
    printf("soudnspd %e, sum_dem %e, P %e, R %e\n",soundspd,sum_dem,C->Flow_data.State_vec.Pres,C->Flow_data.State_vec.Rhof);
    if(C->Merge_data == NULL)
      printf("Ain't merged\n");
    printf("r_grad %e %e %e, E_grad %e %e %e, u_grad %e %e %e, v_grad %e %e %e\n",C->Flow_data.Grads.grad_rhof[0],C->Flow_data.Grads.grad_rhof[1], C->Flow_data.Grads.grad_rhof[2],
	   C->Flow_data.Grads.grad_E[0],C->Flow_data.Grads.grad_E[1], C->Flow_data.Grads.grad_E[2],
	   C->Flow_data.Grads.grad_u[0],C->Flow_data.Grads.grad_u[1],C->Flow_data.Grads.grad_u[2],C->Flow_data.Grads.grad_v[0],
	   C->Flow_data.Grads.grad_v[1],C->Flow_data.Grads.grad_v[2]);
    if(C->Flow_data.Grads.grad_rhof[0] >= 0) {
      if(C->Flow_data.Grads.grad_rhof[0] > DBL_MAX)
	printf("Very +ve\n");
    }
    else if(C->Flow_data.Grads.grad_rhof[0] < 0) {
      printf("A -ve one\n");
      if(C->Flow_data.Grads.grad_rhof[0] < -DBL_MAX)
	printf("Very -ve\n");
    }
    else {
      printf("A truly stuffed drdx\n");
    }
    abort();
  }
  else if (Vol/sum_dem > 0) {
  }
  else {
    printf("Got %g time, C %hd (%10.9e %10.9e %10.9e), type %hd\n",Vol/sum_dem,C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],
	   C->cell_type);
    abort();
  }
#endif

  return(Vol/sum_dem);  
}

/*------------------------------------------------------------------*/

/** \brief Compute exp(-x**2) and erf(x) with a polynomial approximation.
 * \param sn   : IN  : x
 * \param *ex  : OUT : exp(x**2)
 * \param *ef  : OUT : erf(x)  error function
 */
int exxef(double sn, double *ex, double *ef)
{
    double snsq, exx, ef1, y;

#   define P      0.327591100
#   define A1     0.254829592
#   define A2    -0.284496736
#   define A3     1.421413741
#   define A4    -1.453152027
#   define A5     1.061405429
#   define LIMIT  5.0
#   define EXLIM  0.138879e-10
#   define EFLIM  1.0

#   define DSIGN(val,sgn) ( (sgn >= 0.0)? ABS(val): -ABS(val) )

    if (ABS(sn) > LIMIT) {
        exx = EXLIM;
        ef1 = EFLIM;
    } else {
        snsq = sn * sn;
        exx = exp(-snsq);
        y = 1.0 / (1.0 + P * ABS(sn));
        ef1 = 1.0 - y * (A1 + y * (A2 + y * (A3 + y * (A4 + A5 * y)))) * exx;
    }

    *ex = exx;
    *ef = DSIGN(ef1, sn);

    return (0);
}   

/*------------------------------------------------------------------*/

#if SUPERSONIC_VORTEX

double sv_dens(double rhoi, double Mi, double ri, double r)
{
  return(rhoi*pow((1.0 + ((gas_data.gam_amb-1.0)/2.0)*SQR(Mi)*(1.0 - SQR(ri/r))), 1.0/(gas_data.gam_amb-1.0)));
}

short int calc_vortex_err(double rhoi, double Mi, double ri, int run_timesteps_so_far, int num_cells, double current_time)
{
  List_leaf L = Leaves;
  Cart_cell C;
  double max_drr, drr, L1_sum, L2_sum, r, rho_exact;

  max_drr = DBL_MAX;

  while(L != NULL) {
    C = L -> cell_loc;
    
    drr = ABS((C->Flow_data.Temp_state_vec.Rhof) - (C->Flow_data.State_vec.Rhof))/(C->Flow_data.State_vec.Rhof);
    if(drr < max_drr)
      max_drr = drr;

    L = L -> next;
  }

  if(max_drr <= max_drr_tol) { /*We can finish the simulation*/
    printf("\n");
    printf("Now terminating the simulation - max_drr is %g\n",max_drr);

    L1_sum = 0; L2_sum = 0;
    
    L = Leaves;
    while(L != NULL) {
      C = L -> cell_loc;
      
      r = SQRT(SQR(C->centroid[0]) + SQR(C->centroid[1]));
      rho_exact = sv_dens(rhoi, Mi, ri, r);
      
      L1_sum += ABS((C->Flow_data.State_vec.Rhof) - rho_exact)/rho_exact;
      L2_sum += SQR(((C->Flow_data.State_vec.Rhof) - rho_exact)/rho_exact);

      L = L -> next;
    }

    L1_sum = L1_sum/num_cells;
    L2_sum = SQRT(L2_sum/num_cells);
    
    printf("L1 norm is %e, L2 norm is %e\n",L1_sum,L2_sum);
    
    if(write_soln(run_timesteps_so_far, num_cells, current_time) == ERROR)
      printf("   Couldnt' write flow soln at last time 4 some reason\n");

    return(TRUE);
  }

  return(FALSE);
}

#endif
