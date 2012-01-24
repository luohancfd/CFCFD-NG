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


/**\file Source file on adaptation routines for octree cartesian cells - functions associated with cell creation/deletion.
Functions dealing with establishing cell and flux connectivities are dealt with elsewhere.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ov_kernel.h"
#include "ov_adapt.h"
#include "ov_connect.h"
#include "ov_lists.h"
#include "ov_setgeom.h"
#include "ov_icbc.h"
#include "ov_merge.h"
#include "ov_recon.h"
#include "ov_thermo.h"

#define NOTINCLUDE 0
#define DEBUG 1
#define DOPRINTFS 0
#define PARREFS 0

#if !NOTINCLUDE
extern int atime;
Flux_vector sflx = NULL;
extern int stepnow;
#endif

extern short int twoD;
extern Gas_dat gas_data;
extern short int adapt_type;
extern double adapt_param;
extern double adapt_noise;
extern double refine_param;
extern double coarsen_param;
extern short int min_refinement_level;
extern double root_scale;
extern short int max_refinement_level;
extern short int IC_level;
extern double root_shift[3];
extern List_leaf Leaves;
extern List_leaf Leaves_tail;
extern List_leaf * Leaf_heads;
extern List_leaf * Leaf_tails;
extern List_vtx_glob Vtxlist;
extern List_vtx_glob Vtxlist_tail;
extern List_vtx_glob * Vtx_heads;
extern List_vtx_glob * Vtx_tails;
extern short int refining_root_now;
extern short int want_to_adapt;
extern short int adapt_in_parallel;
extern int number_of_threads;

/*------------------------------------------------------------------*/

/**\brief Compute error indicators and flag cells for refinement/coarsening (first step).  Use either Jacob's criterion
   for detecting shocks (velocity differences), or Petrie's density derivatives criterion.*/

void compute_adapt_param(List_leaf L, List_leaf Tail)
{
  List_leaf node = L;
  Cart_cell C;
  short int i, j, k, num_neighbs;
  double shock_param, a, amin, dx;
  double avg_ctr, diff2, diff1, noise, err_param;
  short int use_dir[6] = {FALSE};
  double ctr_rho[6];
  double usgn, vsgn, wsgn;
  num_neighbs = 0; shock_param = 0; amin = 0; dx = 0; avg_ctr = 0;
#if 1
    while(node != NULL) {

      C = node -> cell_loc;

#if DETON
      if(C -> un_det != TRUE)
#endif
	{
	  if(adapt_type == 1) {
	    if(C -> cell_level > min_refinement_level)
	      C -> adapt_flag = COARSEN; /*Let it be coarsenable until proven otherwise for type 1 adaptation type*/
	  }

	  if((adapt_type == 1) || (adapt_type == 3)) {
	    dx = CALC_CELL_EDGE_LENGTH(C);
	    if(gas_data.num_species > 1)
	      amin = get_mixture_soundspd(&(C -> Flow_data.State_vec)); 
	    else amin = SQRT((gas_data.gam_amb)*(C->Flow_data.State_vec.Pres)/(C->Flow_data.State_vec.Rhof));
	  }

	  if((adapt_type == 2) || (adapt_type == 3)) {
	    for(i = 0; i < NUM_FACES; i++) {  /*Reset*/
	      use_dir[i] = FALSE;
	      ctr_rho[i] = 0;
	    }
	  }

	  for(i = 0; i < NUM_FACES; i++) {
	    if((C -> face_neighbours[i][0] != NULL) && (C -> face_neighbours[i][1] == NULL) &&
#if DETON
	       (C -> face_neighbours[i][0] -> un_det != TRUE) &&
#endif
	       ((C->flux_area[i][0] > 0) || (C->flux_area[i][1] > 0) || (C->flux_area[i][2] > 0) || (C->flux_area[i][3] > 0)) &&
	       ((C -> Merge_data == NULL) || (C -> Merge_data != C -> face_neighbours[i][0] -> Merge_data))) { /*1 cell in this direction*/

	      if((adapt_type == 1) || (adapt_type == 3)) {
		if(gas_data.num_species > 1)
		  a = get_mixture_soundspd(&(C -> face_neighbours[i][0] -> Flow_data.State_vec));
		else a = SQRT((gas_data.gam_amb)*(C->face_neighbours[i][0]->Flow_data.State_vec.Pres)/
			      (C->face_neighbours[i][0]->Flow_data.State_vec.Rhof));

		if(a < amin)
		  amin = a;
	      }

	      if((adapt_type == 2) || (adapt_type == 3)) {
		ctr_rho[i] = C -> face_neighbours[i][0] -> Flow_data.State_vec.Rhof;
		use_dir[i] = TRUE;  
	      }
	    }
	    else if(C -> face_neighbours[i][1] != NULL) {
	      if((adapt_type == 2) || (adapt_type == 3)) {
		avg_ctr = 0;
		num_neighbs = 0;
	      }

	      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) {
		if((C -> flux_area[i][j] > 0) && 
#if DETON
		   (C -> face_neighbours[i][j] -> un_det != TRUE) &&
#endif
		   ((C -> Merge_data == NULL) || (C -> Merge_data != C -> face_neighbours[i][j] -> Merge_data))) { /*4 cells here*/

		  if((adapt_type == 1) || (adapt_type == 3)) {
		    if(gas_data.num_species > 1)
		      a = get_mixture_soundspd(&(C -> face_neighbours[i][j] -> Flow_data.State_vec));
		    else a = SQRT((gas_data.gam_amb)*(C->face_neighbours[i][j]->Flow_data.State_vec.Pres)/
				  (C->face_neighbours[i][j]->Flow_data.State_vec.Rhof));
		    if(a < amin)
		      amin = a;
		  }

		  if((adapt_type == 2) || (adapt_type == 3)) {  
		    avg_ctr += C -> face_neighbours[i][j] -> Flow_data.State_vec.Rhof;
		    num_neighbs++;
		  }
		}
	      }

	      if((adapt_type == 2) || (adapt_type == 3)) {
		if(num_neighbs > 0) {
		  ctr_rho[i] = avg_ctr/num_neighbs;
		  use_dir[i] = TRUE;
		}
	      }
	    }
	  }
	  
	  if((adapt_type == 2) || (adapt_type == 3)) {
	    /*Now get numerator and denominator terms for the error indicator*/
	    
	    diff1 = 0;
	    diff2 = 0;
	    noise = 0;

	    for(i = 0; i < 3; i++) {
	      j = 2*i; /*'Positive' directions (E, N, U)*/
	      k = j+1;

	      if((use_dir[j] == TRUE) && (use_dir[k] == TRUE)) { /*Can use this axis*/
		diff2 += ABS(ctr_rho[j] - 2.0*(C -> Flow_data.State_vec.Rhof) + ctr_rho[k]);
		diff1 += ABS(ctr_rho[j] - (C -> Flow_data.State_vec.Rhof)) + ABS((C -> Flow_data.State_vec.Rhof) - ctr_rho[k]);
		noise += ctr_rho[j] + 2.0*(C -> Flow_data.State_vec.Rhof) + ctr_rho[k];
	      }
	    }	  

	    err_param = diff2/(diff1 + (adapt_noise*noise)); /*Error indicator is computed*/

	    if((C -> cell_level < max_refinement_level) && (err_param >= refine_param)) 
	      C -> adapt_flag = REFINE;
	    if((C -> cell_level > min_refinement_level) && (err_param < coarsen_param))
	      C -> adapt_flag = COARSEN;	  	    
	  }

	  if((adapt_type == 1) || ((adapt_type == 3) && (C -> adapt_flag != REFINE))) {
	    for(i = 0; i < 3; i++) {
	      switch(i) {
	      case 0: /*X direction*/
		shock_param = (C->Flow_data.Grads.grad_u[0])*dx/amin; 
		break;
	      case 1:
		shock_param = (C->Flow_data.Grads.grad_v[1])*dx/amin;
		break;
	      case 2:
		shock_param = (C->Flow_data.Grads.grad_w[2])*dx/amin;
		break;
	      }

	      if(shock_param <= -adapt_param) { /*Cells at max level of refinement may not need refinement but also not coarsening*/ 
		C -> adapt_flag = REFINE; /*If adapt_type is 3, this will cancel any previous COARSEN flags*/
		break;
	      }
	    }
	  }

	
	  /*'Protective layer' to ensure shock never passes to coarser cell by refining cells beside it*/
	  if((C -> adapt_flag == REFINE) && (C -> cell_level == max_refinement_level)) {

	    usgn = SGN(C -> Flow_data.State_vec.RhoU);
	    vsgn = SGN(C -> Flow_data.State_vec.RhoV);
	    wsgn = SGN(C -> Flow_data.State_vec.RhoW);
	   	    
	    for(i = 0; i < NUM_FACES; i++) {
	      if(((usgn > 0) && (i == EST)) || ((vsgn > 0) && (i == NTH)) || ((wsgn > 0) && (i == UPR)) ||
		 ((usgn < 0) && (i == WST)) || ((vsgn < 0) && (i == STH)) || ((wsgn < 0) && (i == LWR))) { /*Refine in direction
													     of shock only*/
		if((C -> face_neighbours[i][0] != NULL) && (C -> face_neighbours[i][1] == NULL) &&
#if DETON
		   (C -> face_neighbours[i][0] -> un_det != TRUE) &&
#endif
		   ((C->flux_area[i][0] > 0) || (C->flux_area[i][1] > 0) || (C->flux_area[i][2] > 0) || (C->flux_area[i][3] > 0))) {
		  C -> face_neighbours[i][0] -> adapt_flag = REFINE;
		}
		else if(C -> face_neighbours[i][1] != NULL) {
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) {
		    if((C -> flux_area[i][j] > 0) && 
#if DETON
		       (C -> face_neighbours[i][j] -> un_det != TRUE) 
#endif
		       ) {
		      C -> face_neighbours[i][j] -> adapt_flag = REFINE;
		    }
		  }
		}
	      }
	    }
	  }
	}

      if(node == Tail)
	break;
      else node = node -> next;
    }
#endif
}

/*------------------------------------------------------------------*/

/**\brief Go through list of leaf cells and adapt (in serial/parallel)*/

short int adapt(int *num_refine, int *num_coarsen)
{ 
  short int i;
  Cart_cell C, Parent;
  List_leaf L; 
  *num_refine = 0; *num_coarsen = 0;
  
  if(adapt_in_parallel == FALSE)
    {           
#if NOTINCLUDE
      if(stepnow >= 143)
	printf("About to compute adapt param\n");
#endif

      /*Now compute adaptation parameters for all cells*/
      compute_adapt_param(Leaves, Leaves_tail);
     
      /*Now go through and refine if required*/

#if NOTINCLUDE
      if(stepnow >= 143)
	printf("About to refine\n");
#endif

      L = Leaves;
      while(L != NULL)
	{
	  C = L -> cell_loc;

	  if(C -> adapt_flag == REFINE)
	    {
#if NOTINCLUDE
      if(stepnow >= 145)
	printf("About to refine cell level %hd @ %g %g %g\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
#endif

	      i = refine(C, 1, FALSE, FALSE, &L); /*Refined cells' children put at front of the list*/

#if NOTINCLUDE
      if(stepnow >= 145)
	printf("Finished refinin' cell level %hd @ %g %g %g\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
#endif
	     
	      if(i == ERROR)
		{
		  printf("adapt(): Error refining cell\n");
		  return(ERROR);
		}
	      else if(i == NOERROR)
		(*num_refine)++;
	    } 
	  else L = L -> next;
	}

#if NOTINCLUDE
      if(stepnow >= 143)
	printf("About to coars\n");
#endif

      /*Now go through and coarsen if required - start again from the beginning of list*/
      
      L = Leaves;
      while(L != NULL)
	{ 
	  C = L -> cell_loc; 

	  Parent = C -> parent;

	  while((L != NULL) && (L -> cell_loc -> parent == Parent)) 
	    L = L -> next;

	  if((C -> adapt_flag == COARSEN) || (C -> adapt_flag == REALLY_COARSEN)) 
	    { /*Can coarsen - check other siblings too*/

	      for(i = 0; i < MAX_NUM_CHILDREN; i++) /*See if any siblings can't be coarsened*/
		{
		  if((Parent -> children[i] -> cell_type != SOLID) && 
		     ((Parent -> children[i] -> adapt_flag == UNKNOWN) || (Parent -> children[i] -> adapt_flag == REFINE)))
		    break;
		}

	      if(i == MAX_NUM_CHILDREN) /*So all sibling cells can be coarsened; thus parent can be*/
		{
		  i = coarsen(Parent, C);
			                                                                                  
		  /*Regardless of whether the coarsen operation was successful, we need not visit its
		    children anymore as they've been checked*/
                                                                                
		  if(i == ERROR)
		    {
		      printf("adapt(): Error coarsening cell\n");
		      return(ERROR);
		    }
		  else if(i == NOERROR)
		    (*num_coarsen)++;
		}
	      else  /*Can't coarsen*/
		{
		  for(i = 0; i < MAX_NUM_CHILDREN; i++)
		    Parent -> children[i] -> adapt_flag = UNKNOWN; /*Reset flags anyway*/
		}		      
	    }
	  else /*Not flagged for coarsening*/
	    {
	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Parent -> children[i] -> adapt_flag = UNKNOWN;
	    }
	}

#if NOTINCLUDE
      if(stepnow >= 143)
	printf("About to get new grads and limiters\n");
#endif

      /*Compute gradient and limiter values for newly created cells here*/
      get_grads_all_cells(Leaves, Leaves_tail, TRUE);
      compute_limiter(Leaves, Leaves_tail, TRUE);
      
      /*Allocate flux vectors for new cells*/

#if NOTINCLUDE
      if(stepnow >= 143)
	printf("About to alloc fluxes\n");
#endif

      if(alloc_fluxes(Leaves, Leaves_tail, ALL) == ERROR)
	{
	  printf("adapt(): Error allocating flux vectors\n");
	  return(ERROR);
	}
    }
  else /*Adapt in parallel*/
    { 
      int thread_num = 0;
      int num_ref = 0; 
      int num_coars = 0;      
      int ref, coars;
      short int error_found = FALSE;
      List_leaf *Head; List_leaf *Tail;
      List_vtx_glob *Vhead; List_vtx_glob *Vtail;
      Head = &Leaves;
      Tail = &Leaves_tail;
      Vhead = &Vtxlist;
      Vtail = &Vtxlist_tail;
      coars = 0;
#if 0
      List_leaf L; Cart_cell C; int count, count0, count2, count3, count4, count5, count6, count7;
#endif
#if 0
  L = Leaves;
  count = 0; count0 = 0;
  while(L != NULL) {
    C=L->cell_loc;
    count0++;
    L=L->next;
  }

  L = Leaf_heads[0];
  count = 0; count2=0;
  while(L != NULL) {
    C=L->cell_loc; 
    count2++;
    if(L==Leaf_tails[0])
      break;
    L=L->next;
  }
  
  L = Leaf_heads[1];
  count = 0; count3 = 0;
  while(L != NULL) {
    C=L->cell_loc;
    count3++;
    if(L==Leaf_tails[1])
      break;
    L=L->next;
  }
  count5 = count2+count3;
  printf("\n");
  printf("B4 ad, thus TRUE no. %d, (L1 %d, L2 %d), no. here %d, this&true %d\n",count0,count2,count3,count5,count5==count0);
#endif

      /*Count vtxs here?*/
      count_verticies(); /*Assign verticies to different threads*/

      #pragma omp parallel private(Head, Tail, Vhead, Vtail, thread_num, ref, coars) reduction (+: num_ref, num_coars) 
      {		
        #ifdef _OPENMP
	thread_num = omp_get_thread_num();
        #endif

	Head = &Leaf_heads[thread_num];
	Tail = &Leaf_tails[thread_num]; 
	Vhead = &Vtx_heads[thread_num];
	Vtail = &Vtx_tails[thread_num]; 
	
	compute_adapt_param(*Head, *Tail); /*Go through cells and initially flag their adaptation status*/
        #pragma omp barrier
	
	/*Go through cells and flag them for refinement*/
	flag_for_refine(*Head, *Tail);
        #pragma omp barrier
	
	flag_for_coarsen(*Head, *Tail);
        #pragma omp barrier

#if EXHAUSTIVE_DEBUG
	printout_adapt_cells(REFINE, 0, *Head, *Tail, thread_num, stepnow);
#pragma omp barrier
#endif

        #pragma omp single
	{
	  /*Go through cells in serial and unmerge any that are flagged for refinement/coarsening*/
	  unmerge_all_cells(); 

	  modify_lists(Leaf_heads, Leaf_tails); /*Ensure no siblings to be coarsened are split up between different lists*/
	  break_list(Leaf_heads, Leaf_tails); /*Break up global list here*/
	  break_vtx_list(Vtx_heads, Vtx_tails);
	}
        #pragma omp barrier

#if 0
  L = Leaf_heads[thread_num];
  count = 0;
  while(L != NULL) {
    C=L->cell_loc;
    if((C->cell_level==7)&&(C->centroid[0]<=0.015625)&&(C->centroid[1]>=-0.015625)&&(C->centroid[2]>=-0.171875)&&
       (C->centroid[2]<=-0.15625)) {
      printf("After modifyin t %d - We've found child %hd (%e %e %e)\n",thread_num,C->child_num,C->centroid[0],C->centroid[1],C->centroid[2]);
      count++;
    }
      
      L=L->next;
    }
  if(count == 0)
    printf("t %d, After modifying found nothin'\n",thread_num);
#pragma omp barrier
#endif

	/*Begin housekeeping steps for cell refinement (5 steps)*/	
	ref = par_refine(Head, Tail, Vhead, Vtail, ALLOC_CHILD);
	if(ref == ERROR) /*Or can use NULL*/
	  {
	    error_found = TRUE;
	    #pragma omp flush(error_found)
	  }
	num_ref = num_ref + ref;
        #pragma omp barrier

	if(error_found == FALSE)
	  if(par_refine(Head, Tail, Vhead, Vtail, POSITIVE) == ERROR)
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }
        #pragma omp barrier

	if(error_found == FALSE)
	  if(par_refine(Head, Tail, Vhead, Vtail, NEGATIVE) == ERROR)
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }
        #pragma omp barrier

	if(error_found == FALSE)
	  if(par_refine(Head, Tail, Vhead, Vtail, POSITIVE2) == ERROR)
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }
        #pragma omp barrier

	if(error_found == FALSE)
	  if(par_refine(Head, Tail, Vhead, Vtail, NEGATIVE2) == ERROR)
	    {
	      error_found = TRUE;
              #pragma omp flush(error_found)
	    }
        #pragma omp barrier

#if EXHAUSTIVE_DEBUG
	printout_adapt_cells(REALLY_COARSEN, 0, *Head, *Tail, thread_num, stepnow);
#pragma omp barrier
#endif

	/*Do necessary housekeeping for soon-to-be coarsened cells (in 4 steps)*/

	if(error_found == FALSE)
	  coars = par_coarsen(Head, Tail, Vhead, Vtail, POSITIVE);
	  if(coars == ERROR)
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }
	num_coars = num_coars + coars;
        #pragma omp barrier

	if(error_found == FALSE)
	  if(par_coarsen(Head, Tail, Vhead, Vtail, NEGATIVE) == ERROR)
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }
        #pragma omp barrier

	if(error_found == FALSE)
	  if(par_coarsen(Head, Tail, Vhead, Vtail, POSITIVE2) == ERROR)
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }
        #pragma omp barrier

	if(error_found == FALSE)
	  if(par_coarsen(Head, Tail, Vhead, Vtail, NEGATIVE2) == ERROR)
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }
        #pragma omp barrier
	
	/*Delete verticies that need to be deleted*/

	if(error_found == FALSE)
	  delete_verticies(Vhead, Vtail); /*Don't need to put a barrier here*/ 

	/*Compute gradient and limiter values for newly created cells here*/

	if(error_found == FALSE) 
	  {
	    get_grads_all_cells(*Head, *Tail, TRUE);
	    compute_limiter(*Head, *Tail, TRUE);
	  }

	/*Allocate fluxes*/

	if(error_found == FALSE)
	  if(alloc_fluxes(*Head, *Tail, POSITIVE) == ERROR) /*Can use NULL instead of Leaves_tail*/
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }

        #pragma omp barrier

	if(error_found == FALSE)
	  if(alloc_fluxes(*Head, *Tail, NEGATIVE) == ERROR)
	    {
	      error_found = TRUE;
	      #pragma omp flush(error_found)
	    }

#if EXHAUSTIVE_DEBUG
#pragma omp barrier
	printout_adapt_cells(1, 1, *Head, *Tail, thread_num, stepnow);
#endif
      }
      
      /*Join up sublists again*/
    
 
      merge_lists(Leaf_heads, Leaf_tails);
      merge_vtx_lists(Vtx_heads, Vtx_tails);
      Leaves = Leaf_heads[0];
      Leaves_tail = Leaf_tails[number_of_threads-1];
      Vtxlist = Vtx_heads[0];
      Vtxlist_tail = Vtx_tails[number_of_threads-1];
	
      *num_refine = num_ref;
      *num_coarsen = num_coars;

#if EXHAUSTIVE_DEBUG
L = Leaves;
while(L != NULL) {
C = L->cell_loc;
C->just_created = FALSE;

L = L -> next;
}
#endif

      if(error_found == TRUE)
	return(ERROR);
    }  

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Go through list (or sublist) and flag a cell for refinement*/

void flag_for_refine(List_leaf L, List_leaf Tail)
{
  List_leaf node = L;
  Cart_cell C;

  while(node != NULL)
    {
      C = node -> cell_loc;

      if(C -> adapt_flag == REFINE)
	{
	  if(C -> cell_level < max_refinement_level) /*Must test if maximum level not exceeded*/
	    check_neighbour_for_refine(C, NULL);
	  else C -> adapt_flag = UNKNOWN; /*Can't refine*/
	}

      if(node == Tail)
	break;
      else node = node -> next;
    }
}

/*------------------------------------------------------------------*/

/**\brief Go through list (or sublist) and flag a cell for coarsening*/

void flag_for_coarsen(List_leaf L, List_leaf Tail)
{
  short int i;
  List_leaf node = L;
  Cart_cell C, Parent;

  while(node != NULL)
    {
      C = node -> cell_loc;
      Parent = C -> parent;

      if(C -> adapt_flag == COARSEN)
	{
	  for(i = 0; i < MAX_NUM_CHILDREN; i++) 
	    {
	      if((Parent -> children[i] -> cell_type != SOLID) && (Parent -> children[i] -> adapt_flag != COARSEN) &&
		 (Parent -> children[i] -> adapt_flag != REALLY_COARSEN))
		break;
	    }

	  if(i == MAX_NUM_CHILDREN)
	    check_cell_for_coarsening(Parent); /*Will check if parent can really be coarsened, then will set 
						 adapt_flags of children to REALLY_COARSEN*/
	  else
	    {
	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		{
		  if(Parent -> children[i] -> adapt_flag != REFINE)
		    Parent -> children[i] -> adapt_flag = UNKNOWN;
		} 
	    }
	}

      if(node == Tail)
	break;
      else node = node -> next;
    }
}

/*------------------------------------------------------------------*/

/**\brief Refine cell in parallel (done in 5 stages)*/

int par_refine(List_leaf *L, List_leaf *Tail, List_vtx_glob *Vtxhead, List_vtx_glob *Vtxtail, short int stage)
{
  List_leaf node = *L;
  Cart_cell C;
  short int i;
  int num_refine = 0;

  short int rm_tmp_list = FALSE;

  if(stage == NEGATIVE)
    rm_tmp_list = TRUE;

  while(node != NULL)
    {
      C = node -> cell_loc; /*Traverse list and perform housekeeping for any cells that need to be refined*/

      if(C -> adapt_flag == REFINE)
	{ 
	  if(stage == ALLOC_CHILD) /*1st stage - allocate children, do housekeeping that's independent of neighbours*/
	    {
	      if(make_children(C) == ERROR)
		return(ERROR);
	      
	      if(assign_vtx_after_refine(C -> children, stage, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      if(set_cell_geom_data(C -> children, 1, TRUE, TRUE, FALSE, ALL, FALSE) == ERROR)
		return(ERROR);
		
	      reconstruct_refined_children(C -> children);

	      num_refine++; /*Increase counters*/
	
#if 0
	      printf("We're refining cell %hd (%e %e %e)\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
#endif

	      node = node -> next;
	    }	
	  else if((stage == POSITIVE) || (stage == NEGATIVE))
	    { 			
	      establish_neighbour_after_refine(C -> children, stage);

	      if(assign_vtx_after_refine(C -> children, stage, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);
	      
	      if(set_cell_geom_data(C -> children, 1, FALSE, FALSE, TRUE, stage, rm_tmp_list) == ERROR)
		return(ERROR);
	  	  
	      reassign_fluxes_after_refine(C -> children, stage);	      
	      
	      node = node -> next;		
	    }
	  else if(stage == POSITIVE2)
	    {
	      if(assign_vtx_after_refine(C->children, stage, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      node = node -> next;
	    }
	  else if(stage == NEGATIVE2)
	    {
	      if(assign_vtx_after_refine(C -> children, stage, Vtxhead, Vtxtail) == ERROR) /*Still do other phases of vertex assigning*/
		return(ERROR);

	      /*Can finally add children to list and delete parent (at this stage it won't affect other 
		housekeeping operations)*/

	      C -> adapt_flag = UNKNOWN; /*Reset at the end*/

	      node = node -> next; /*Move to next pointer before deleting*/
			
	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		{
#if NOTINCLUDE /*Check if this cell is the only cell by itself, if so, then it's irrelevant - but this isn't a very 
		 pratical possibility*/ 
		  j = 0; k = 0;
		  
		  if(C -> children[i] -> cell_type != SOLID)
		    {
		      for(j = 0; j < NUM_FACES; j++) 
			{
			  if((C -> children[i] -> face_neighbours[j][0] != NULL) && (C -> children[i] -> face_neighbours[j][1] == NULL) && 
			     ((C -> children[i] -> flux_area[j][0] != 0) || (C -> children[i] -> flux_area[j][1] != 0) || 
			      (C -> children[i] -> flux_area[j][2] != 0) || (C -> children[i] -> flux_area[j][3] != 0)))
			    break;
			  else if(C -> children[i] -> face_neighbours[j][1] != NULL) 
			    {
			      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) 
				if(C -> children[i] -> flux_area[j][k] != 0)  
				  break;

			      if(k != MAX_NUM_NEIGHBOURS)
				break;			    
			    }
			}
		    }

		  if((C -> children[i] -> cell_type == SOLID) || ((j*k) == (NUM_FACES*MAX_NUM_NEIGHBOURS)))
		    ; /*Don't add to leaf list - Can't find any neighbour that can flux/merge with - cell's solid/by itself (irrelevant)*/
#endif

		  if(C -> children[i] -> cell_type != SOLID)
		    if(add_to_leaf_list_front(C -> children[i], &(C -> Leaf_list_loc), L, Tail) == ERROR) 
		      return(ERROR);		    
		}

	      if(C -> Leaf_list_loc != NULL)
		delete_from_leaf_list(C, L, Tail);
	    }	    
	}
      else node = node -> next;
    }

  if(stage == ALLOC_CHILD)
    return(num_refine);
  else return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Coarsen cells in parallel (done in 4 steps)*/

int par_coarsen(List_leaf *L, List_leaf *Tail, List_vtx_glob *Vtxheads, List_vtx_glob *Vtxtails, short int stage)
{
  List_leaf node = *L;
  Cart_cell C, Parent;
  short int i, j, k, x;
  int num_coarsen = 0;
  k = 0;
  while(node != NULL)
    { 
      C = node -> cell_loc;
      Parent = C -> parent;
	  
      if(C -> adapt_flag == REALLY_COARSEN)
	{ 
	  if((stage == POSITIVE) || (stage == NEGATIVE))
	    {
	      establish_neighbour_after_coarsen(Parent, stage);

	      if(stage == POSITIVE)
		{
		  reconstruct_coarsened_cell(Parent);
		  
		  num_coarsen++;
		}

	      reassign_fluxes_after_coarsen(Parent -> children, stage);

	      assign_vtx_after_coarsen(Parent -> children, stage, Vtxheads, Vtxtails);
		
	      while((node != NULL) && (node -> cell_loc -> parent == Parent)) 
		node = node -> next;
	    }
	  else if(stage == POSITIVE2)
	    {
	      assign_vtx_after_coarsen(Parent -> children, stage, Vtxheads, Vtxtails);
	      
	      while((node != NULL) && (node -> cell_loc -> parent == Parent))
		node = node -> next;	      
	    }
	  else if(stage == NEGATIVE2)
	    {	
	      assign_vtx_after_coarsen(Parent -> children, stage, Vtxheads, Vtxtails);
	      	      
	      /*We must update the areas again due to solid degeneracies where a newly solid neighbour cell updates with
		a parent's children (but check_cell_for_coarsening() is done before par_refine())*/
	      
	      for(i = 0; i < NUM_FACES; i++) /*Cycle through each face*/
		{
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Cycle through each quadrant*/
		    {
		      Parent -> flux_area[i][j] = 0;
                                                                                                                           
		      switch(i)
			{
			case NTH:
			  k = 2*j+1; break; /*kth child cell for this quadrant*/
			case STH:
			  k = 2*j; break;
			case EST:
			  if(j <= 1)
			    k = j+2;
			  else k = j+4;
			  break;
			case WST:
			  if(j <= 1)
			    k = j;
			  else k = j+2; break;
			case LWR:
			  k = j; break;
			case UPR:
			  k = j+4; break;
			}
                                                                                                                           
		      for(x = 0; x < MAX_NUM_NEIGHBOURS; x++)
			Parent -> flux_area[i][j] += (Parent -> children[k] -> flux_area[i][x]);
		    }
		}

#if 0
	      printf("We're coarsen cell %hd (%e %e %e)\n",Parent->cell_level,Parent->centroid[0],Parent->centroid[1],Parent->centroid[2]);
#endif

		/*Now can remove children and update lists - at this stage all coarsening housekeeping 
		  operations will be thread-safe */

	      if(Parent -> Leaf_list_loc == NULL)
		if(add_to_leaf_list_front(Parent, &node, L, Tail) == ERROR)
		  return(ERROR);

	      while((node != NULL) && (node -> cell_loc -> parent == Parent)) /*Skip ahead to next set of sibling cells*/
		node = node -> next;

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		{
		  if(Parent -> children[i] -> Leaf_list_loc != NULL)
		    delete_from_leaf_list(Parent -> children[i], L, Tail);
		  
		  free(Parent -> children[i] -> cell_num);
		  free(Parent -> children[i]);
		  Parent -> children[i] = NULL;
		}

	      Parent -> just_created = TRUE;
	    }      
	}
      else 
	{
	  while((node != NULL) && (node -> cell_loc -> parent == Parent))
	    node = node -> next;
	}
    }  
    
  if(stage == POSITIVE)
    return(num_coarsen);
  else return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Refine cell (and its neighbours if required - to ensure grid smoothness) num_times in serial.  Parallel refinement
 requires somewhat different treatment*/

short int refine(Cart_cell C, short int num_times, short int additionally_refine, short int set_IC, List_leaf * current_node) 
{ /**< If additionally_refine TRUE, will refine descendants if needed, if set_IC TRUE, will set IC_flags of children*/

  short int i;

  if((C -> children[0]) == NULL) /*Children not already allocated*/
    {	  
      if((C -> cell_level < max_refinement_level) || (refining_root_now == TRUE))
	{ 	     
	  /*Check neighbours if any need refinement (to ensure grid smoothness)*/
      
	  Cart_cell refine_array[3]; /*At most 3 neighbours should be refined*/
	  short int n;
      
	  n = check_neighbour_for_refine(C, refine_array); 

	  if(n > 0) /*Refine neighbours first*/
	    { 
	      for(i = 0; i < n; i++)
		{
		  if(refine(refine_array[i], 1, FALSE, set_IC, &(*current_node)) == ERROR) /*Refine neighbours only once*/
		    return(ERROR);  
		}
	    }

	  if(make_children(C) == ERROR)
	    return(ERROR);
 
	  if(((*current_node) != NULL) && ((*current_node) == C -> Leaf_list_loc))
	    *current_node = (*current_node) -> next; 

	  establish_neighbour_after_refine(C -> children, ALL); /*Update neighbour connectivities*/

#if VTXM
	  if(assign_vtx_after_refine(C -> children, ALL, &Vtxlist, &Vtxlist_tail) == ERROR) /*Assign vertex pointers*/
	    return(ERROR);
#endif

	  if(C -> Merge_data != NULL) /*Unmerge the cell*/
	    unmerge_cell(C, FALSE);
	  
	  if(set_cell_geom_data(C -> children, num_times, TRUE, TRUE, TRUE, ALL, TRUE) == ERROR) /*Set geometric properties of cells*/	
	    return(ERROR);
	  	  
	  if(set_IC != TRUE)
	    reassign_fluxes_after_refine(C -> children, ALL); /*Reassign flux interface pointers from parent to children*/

	  if(set_IC != TRUE)
	    reconstruct_refined_children(C -> children); /*Reconstruct flow variables here (have to do so each time)*/
	  else set_IC_flag(C -> children); /*No need to reconstruct during initial refine of grid*/
	  
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if(C -> children[i] -> cell_type != SOLID)  
		if(add_to_leaf_list_front(C -> children[i], &(C -> Leaf_list_loc), &Leaves, &Leaves_tail) == ERROR) 
		  return(ERROR);	      	      
	    }

	  if(C -> Leaf_list_loc != NULL) /*Only now delete the list node*/
	    delete_from_leaf_list(C, &Leaves, &Leaves_tail); 	    

	  if(num_times == 1)
	    {		
	      for(i = 0; i < MAX_NUM_CHILDREN; i++) 
		{ 
		  if((C -> children[i] -> additional_refine_times == 0) && 
		     ((set_IC == FALSE) || 
#if FLOAT_EQ_STABLE
		      EQ(ROUND(C -> children[i] -> IC_flag), (C -> children[i] -> IC_flag))
#else
		      (ROUND(C -> children[i] -> IC_flag) == C -> children[i] -> IC_flag)
#endif 
		      || (C -> children[i] -> cell_level >= IC_level))) /*If IC_flags rounded are themselves, isn't intersected*/
		    { /*Cell doesn't need further refinement - put on list*/
		      
		      if(set_IC == TRUE)
			set_cell_IC(C -> children[i]); 
		    }
		  else /*Cell needs further refinement*/ 
		    { 				  
		      if(additionally_refine == TRUE)
			{
			  /*Refine for combination of IC/geometrical reasons*/	
#if !NOTINCLUDE
			  if((C -> children[i] -> cell_type != SOLID) && 
			     ((C -> children[i] -> additional_refine_times > 0) || 
			      ((set_IC == TRUE) &&
#if FLOAT_EQ_STABLE
			       NOT_EQ(ROUND(C -> children[i] -> IC_flag), (C -> children[i] -> IC_flag))
#else
			       (ROUND(C -> children[i] -> IC_flag) != C -> children[i] -> IC_flag) 
#endif
			       && (C -> children[i] -> cell_level < IC_level)))) {

			    if(refine(C -> children[i], MAX(1, C -> children[i] -> additional_refine_times), 
				      additionally_refine, set_IC, &(*current_node)) == ERROR)
			      return(ERROR);
			  }
#else	  
			  if((set_IC == TRUE) && (ROUND(C -> children[i] -> IC_flag) != C -> children[i] -> IC_flag) &&  
			     (C -> children[i] -> cell_level < IC_level) && (C -> children[i] -> cell_type != SOLID))
			    {
			      if(refine(C -> children[i], MAX(1, C -> children[i] -> additional_refine_times), 
					additionally_refine, set_IC, &(*current_node)) == ERROR)
				return(ERROR);
			    } 
			  else if((C -> children[i] -> cell_type != SOLID) && (C -> children[i] -> additional_refine_times > 0))
			    {
			      if(refine(C -> children[i], C -> children[i] -> additional_refine_times, 
					additionally_refine, set_IC, &(*current_node)) == ERROR)
				return(ERROR);
			      
			      /*Refine only for geometrical reasons*/
			    }
#endif
			}
		    }
		} 
	    }
	  else if(num_times > 1) /*Haven't reached 'leaf' level of refinement yet*/
	    {
	      num_times--;
	      
	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		{			      
		  if(C -> children[i] -> cell_type != SOLID) 
		    {
		      if(refine(C -> children[i], num_times, additionally_refine, set_IC, &(*current_node)) == ERROR)
			return(ERROR);
		    } /*Refine each child*/
		}
	    } /*Geometrical data successfully set*/
	  
	  C -> adapt_flag = UNKNOWN;
	}
      else
	{
	  /*Move pointers along even if cell can't be refined*/
                                                                                                                             
	  if(((*current_node) != NULL) && ((*current_node) == C -> Leaf_list_loc))
	    *current_node = (*current_node) -> next;
                                   
	  C -> adapt_flag = UNKNOWN; /*Reset any 'REFINE' flags during flow adaptation*/
                                                                                                    
	  return(ERROR2); /*Maximum permissible level exceeded*/
	}     
    }        
  else /*Cell already has children allocated*/
    {	
      if(num_times == 1)
	{
	  /*Note, C -> additional_refine_times has also been computed for already-refined cells (and also IC_flags)*/
	      
	  if(additionally_refine == TRUE) /*Further refinement is necessary*/
	    {
	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		{	
#if !NOTINCLUDE
		  if((C -> children[i] -> cell_type != SOLID) && 
		     ((C -> children[i] -> additional_refine_times > 0) || 
		      ((set_IC == TRUE) && 
#if FLOAT_EQ_STABLE
		       NOT_EQ(ROUND(C -> children[i] -> IC_flag), (C -> children[i] -> IC_flag))
#else
		       (ROUND(C -> children[i] -> IC_flag) != C -> children[i] -> IC_flag) 
#endif
		       && (C -> children[i] -> cell_level < IC_level))))
		    if(refine(C -> children[i], MAX(1, C -> children[i] -> additional_refine_times), 
			      additionally_refine, set_IC, &(*current_node)) == ERROR)
		      return(ERROR);
#else
		  if((set_IC == TRUE) && (ROUND(C -> children[i] -> IC_flag) != C -> children[i] -> IC_flag) && 
		     (C -> children[i] -> cell_level < IC_level) && (C -> children[i] -> cell_type != SOLID)) 
		    {
		      if(refine(C -> children[i], MAX(1, C -> children[i] -> additional_refine_times), 
				additionally_refine, set_IC, &(*current_node)) == ERROR)
			return(ERROR);
		    }
		  else if((C -> children[i] -> cell_type != SOLID) && (C -> children[i] -> additional_refine_times > 0))
		    {
		      if(refine(C -> children[i], C -> children[i] -> additional_refine_times, 
				additionally_refine, set_IC, &(*current_node)) == ERROR)
			return(ERROR); 
		    }
#endif
		}
	    }
	}
      else if(num_times > 1) /*Haven't reached 'leaf' level of refinement yet*/
	{
	  num_times--;
	      
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if(C -> children[i] -> cell_type != SOLID)
		{
		  if(refine(C -> children[i], num_times, additionally_refine, set_IC, &(*current_node)) == ERROR)
		    return(ERROR); 
		  /*Refine each child*/
		}
	    }
	}	     
    }
  
  return(NOERROR); 
}

/*------------------------------------------------------------------*/

/**\brief Coarsen a cell in serial with all non-solid leaf children.  Will check (fluxable) neighbour level differences
   if coarsening allowed.  Parallel coarsening requires somewhat different treatment.*/

short int coarsen(Cart_cell C, Cart_cell child)
{ 
  short int i;
      	  
  i = check_cell_for_coarsening(C);
	  
  if(i == ERROR)
    return(ERROR);
  else if(i == FALSE)
    return(ERROR2); /*Can't coarsen cell*/

  if(C -> Leaf_list_loc == NULL)
    if(add_to_leaf_list_front(C, &(child -> Leaf_list_loc), &Leaves, &Leaves_tail) == ERROR)
      return(ERROR);

  for(i = 0; i < MAX_NUM_CHILDREN; i++)
    {
      if(C -> children[i] -> Leaf_list_loc != NULL) /*They should be on the list anyway*/
	delete_from_leaf_list(C -> children[i], &Leaves, &Leaves_tail);
    }
			  
  establish_neighbour_after_coarsen(C, ALL); /*Now update neighbour connectivities*/
			  
  reconstruct_coarsened_cell(C); /*Before freeing children, reconstruct state vector for cell*/
	    
  for(i = 0; i < MAX_NUM_CHILDREN; i++)
    {
      if(C -> children[i] -> Merge_data != NULL) /*A child is merged - 'unmerge' it*/
	unmerge_cell(C -> children[i], FALSE);
    }

  reassign_fluxes_after_coarsen(C -> children, ALL); /*After successful coarsening, reassign flux vectors*/
      
#if VTXM
  assign_vtx_after_coarsen(C -> children, ALL, &Vtxlist, &Vtxlist_tail); /*Reassign vertex pointers*/
#endif
      
  for(i = 0; i < MAX_NUM_CHILDREN; i++)
    {							  
      free(C -> children[i] -> cell_num);
      free(C -> children[i]);
      C -> children[i] = NULL; /*Finally free the children if required*/
    }
  
  C -> just_created = TRUE;

  return(NOERROR); /*Coarsening operation finally successful*/
}

/*------------------------------------------------------------------*/

/**\brief Allocate space for a single root cell and initialize structure variables appropriately*/

Cart_cell make_root(void) 
{
  Cart_cell root = malloc(sizeof(struct cart_cell)); 

  if((root != NULL) && (add_to_leaf_list_front(root, &Leaves, &Leaves, &Leaves_tail) == NOERROR)) 
    {	
      /*We only want to initialize relevant variables*/

      short int i, j;

      /*No need to Flow_data variables; they're never used anyway*/

      root -> child_num = -1; /*Create special number for root cell*/
      
#if SOLNSPLIT
      root -> cell_num = malloc(sizeof(char));
      if(root -> cell_num == NULL)
	return(NULL);
      else root -> cell_num[0] = '0'; /*Root is uniquely cell_num 0*/
#endif

      for(i = 0; i < 3; i++)
	root -> centroid[i] = root_shift[i];
      
      root -> cell_level = 0; 
      
      root -> cell_length = root_scale;

      root -> parent = NULL; 
      
      for(i = 0; i < MAX_NUM_CHILDREN; i++)
	root -> children[i] = NULL; 
            
      for(i = 0; i < NUM_FACES; i++)
	{
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
	    {
	      root -> face_neighbours[i][j] = NULL; 
	      root -> flux_area[i][j] = UNKNOWN; 
	    }
	}

      /*Assign corner verticies*/

#if VTXM
      for(i = 0; i < MAX_NUM_CHILDREN; i++)
	{
	  root -> verticies[i] = malloc(sizeof(struct vertex));
	  if(root -> verticies[i] == NULL)
	    {
	      for(j = 0; j < i; j++)
		free(root -> verticies[i]);
	      return(NULL);
	    }
	  else
	    {
	      root -> verticies[i] -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(root -> verticies[i] -> vtx_nums == NULL)
		return(NULL);

	      if(add_to_Vtxlist(root -> verticies[i], &Vtxlist, &Vtxlist_tail) == ERROR)
		return(NULL);
	      	      
	      get_vtx_loc(root, i, root -> verticies[i] -> loc);

	      for(j = 0; j < MAX_NUM_CHILDREN; j++)
		root -> verticies[i] -> Leaf_cells[j] = NULL;

	      root -> verticies[i] -> Leaf_cells[7-i] = root;		
	    }
	}
#endif
      
      /*It might be pointless to initialize some variables below, but never mind*/

      
#if AXI_OK
      if(twoD == 2) {
	root -> r_centroid = UNKNOWN;
	for(i = 0; i < 2; i++)
	  for(j = 0; j < 2; j++)
	    root -> face_r_centroid[i][j] = 0;
      }
#endif

      root -> Flow_data.Obstructed_flux = NULL;
      root -> Merge_data = NULL;
      root -> Merge_list_loc = NULL;
      root -> wall_area = UNKNOWN;
      root -> cell_volume = UNKNOWN;
      root -> is_small_cell = UNKNOWN; 
      root -> cell_type = UNKNOWN; 
      root -> IC_flag = UNKNOWN;
      root -> un_det = UNKNOWN; 
      root -> Solid = NULL; 
      root -> additional_refine_times = UNKNOWN; 
      root -> adapt_flag = UNKNOWN;
      root -> shocked_cell = 'n';
      root -> cell_marked = UNKNOWN;
      root -> just_created = TRUE;
      root -> Intersected_bbox = NULL;
      root -> bbox_unknown = TRUE;
    }  
  else 
    {
      if(root != NULL)      
	free(root);	  
      root = NULL;
    }	
  
  return(root);
}

/*------------------------------------------------------------------*/

/**\brief Allocate space for children of a cell and initialize structure variables appropriately*/

short int make_children(Cart_cell C) 
{                               
  short int i, j, k; 
  
  for(i = 0; i < MAX_NUM_CHILDREN; i++)
    {
      C -> children[i] = malloc(sizeof(struct cart_cell));

      if(C -> children[i] == NULL)
	{
	  for(j = 0; j < i; j++)
	    {
	      free(C -> children[j]);
	      C -> children[j] = NULL;
	    }
	  break;
	}      
    }

  /*So all children successfully allocated*/

  if(i == MAX_NUM_CHILDREN) 
    {
      double centroids[MAX_NUM_CHILDREN][3]; /*Centroids of all children*/
      compute_children_centroids(C, centroids);

      for(i = 0; i < MAX_NUM_CHILDREN; i++) 
	{ 
	  for(j = 0; j < NUM_FACES; j++)
	    {
	      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
		{
		  C -> children[i] -> Flow_data.Face_fluxes[j][k] = NULL;
		  C -> children[i] -> face_neighbours[j][k] = NULL;
		} 
	    }
	
	  C -> children[i] -> Flow_data.Obstructed_flux = NULL;
	  
	  if((want_to_adapt == TRUE) || (want_to_adapt == FALSE))
	    {
	      for(j = 0; j < 3; j++)
		C -> children[i] -> centroid[j] = centroids[i][j];

	      C -> children[i] -> child_num = i;
	      
	      C -> children[i] -> cell_level = (C -> cell_level) + 1;

	      C -> children[i] -> cell_length = (C -> cell_length)/2.0;
	      
	      C -> children[i] -> parent = C; 
	      	     
#if SOLNSPLIT
	      C -> children[i] -> cell_num = malloc(((C -> children[i] -> cell_level)+1)*sizeof(char));
	      strncpy(C -> children[i] -> cell_num, "\0", C -> children[i] -> cell_level + 1);
	      if(C -> children[i] -> cell_num == NULL)
		return(ERROR);
	      else {
		for(j = 0; j < ((C -> cell_level) + 1); j++)
		  C -> children[i] -> cell_num[j] = C -> cell_num[j]; /*Only last digit is uniquely the child's*/
		switch(i) {
		case 0:
		  C -> children[i] -> cell_num[(C->cell_level)+1] = '0'; break;
		case 1:
		  C -> children[i] -> cell_num[(C->cell_level)+1] = '1'; break;
		case 2:
		  C -> children[i] -> cell_num[(C->cell_level)+1] = '2'; break;
		case 3:
		  C -> children[i] -> cell_num[(C->cell_level)+1] = '3'; break;
		case 4:
		  C -> children[i] -> cell_num[(C->cell_level)+1] = '4'; break;
		case 5:
		  C -> children[i] -> cell_num[(C->cell_level)+1] = '5'; break;
		case 6:
		  C -> children[i] -> cell_num[(C->cell_level)+1] = '6'; break;
		case 7:
		  C -> children[i] -> cell_num[(C->cell_level)+1] = '7'; break;
		}
	      }
#endif
 
	      for(j = 0; j < MAX_NUM_CHILDREN; j++)
		C -> children[i] -> children[j] = NULL; 
	      
	      C -> children[i] -> cell_volume = UNKNOWN;
	      
	      for(j = 0; j < NUM_FACES; j++)
		{
		  for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
		    C -> children[i] -> flux_area[j][k] = UNKNOWN; 
		}

	      C -> children[i] -> wall_area = UNKNOWN;

	      for(j = 0; j < 3; j++)
		C -> children[i] -> wall_norm[j] = UNKNOWN;
	      
	      C -> children[i] -> cell_type = UNKNOWN; 
	      
	      C -> children[i] -> IC_flag = UNKNOWN;

	      if(C -> un_det == TRUE)
		C -> children[i] -> un_det = TRUE;
	      else C -> children[i] -> un_det = UNKNOWN;
	      
	      C -> children[i] -> additional_refine_times = UNKNOWN;
	      
	      C -> children[i] -> Solid = NULL; 
	      
	      C -> children[i] -> Leaf_list_loc = NULL;

	      C -> children[i] -> is_small_cell = UNKNOWN;

	      C -> children[i] -> Merge_data = NULL;

	      C -> children[i] -> Merge_list_loc = NULL;

	      for(j = 0; j < MAX_NUM_CHILDREN; j++)
		C -> children[i] -> verticies[j] = NULL;

	      /*Pointers to list nodes on vertex lists need not be allocated - they depend on pointers to verticies*/

	      C -> children[i] -> adjacent_to_normal_cell = UNKNOWN;

	      C -> children[i] -> adapt_flag = UNKNOWN;

	      C -> children[i] -> shocked_cell = 'n';

	      C -> children[i] -> cell_marked = UNKNOWN;

	      C -> children[i] -> just_created = TRUE;

	      C -> children[i] -> Intersected_bbox = NULL;

	      C -> children[i] -> bbox_unknown = TRUE;

#if AXI_OK
	      C -> children[i] -> r_centroid = UNKNOWN;
	      for(j = 0; j < 2; j++) {
		C -> children[i] -> face_r_centroid[j][0] = 0;
		C -> children[i] -> face_r_centroid[j][1] = 0;
	      }
#endif

	      /*It might be good to initialize these values anyway*/
	      C -> children[i] -> Flow_data.State_vec.Rhof = 0;
	      C -> children[i] -> Flow_data.State_vec.RhoU = 0;
	      C -> children[i] -> Flow_data.State_vec.RhoV = 0;
	      C -> children[i] -> Flow_data.State_vec.RhoW = 0;
	      C -> children[i] -> Flow_data.State_vec.RhoE = 0;
	      C -> children[i] -> Flow_data.State_vec.rhof_prod = 0;
	      C -> children[i] -> Flow_data.State_vec.Pres = 0;
	      C -> children[i] -> Flow_data.State_vec.T = 0;
#if AFTERBURN
	      C -> children[i] -> Flow_data.State_vec.rhof_prod_ub = 0;
#endif
	    }
	  else 
	    {
	      free(C -> children[i]);
	      C -> children[i] = NULL;
	      return(ERROR);
	    }
	}

      return(NOERROR);       
    } 
  else return(ERROR);
}

/*------------------------------------------------------------------*/
/**\brief Very similar to the serial refine(), except here we want to refine a specific cell
   given by its cell_num.  Useful for reconstructing a mesh for solution continuation, so is a simplier function than refine()*/

short int refine_specific(Cart_cell C, short int cell_level, char *cell_num, Cart_cell * this_cell)
{
  short int i, j;

  if(C -> children[0] != NULL) {
#if 0
    if((cell_level==0) && (cell_num[0]=='0')) {
      printf("We're lookin 4 cell 0 into an already cell with num %s\n",C->cell_num);
    }
#endif

    for(i = 0; i < MAX_NUM_CHILDREN; i++) {
      /*Which part of the tree should we traverse down?  We need to compare all digits of the cell no. up to the 
	nth digit, where n is the child's level*/

      for(j = 0; j < ((C -> children[0] -> cell_level)+1); j++) 
	if(C -> children[i] -> cell_num[j] != cell_num[j])
	  break;
      
      if(j == (C -> children[0] -> cell_level+1)) { /*The cell is down (or is) this branch*/
#if 0
	if((cell_level==0) &&(cell_num[0]=='0')) {
	  printf("We're lookin 4 cell 0 into an already child with num %s\n",C->children[0]->cell_num);
	}
#endif

	refine_specific(C -> children[i], cell_level, cell_num, &(*this_cell));
	break;
      }
    }
  }
  else { /*So we want to refine this cell.  Proceed to do a number of things similar to refine()*/

    /*Notably we don't refine neighbours as all cells will be listed (and sometimes refining neighbours unnecessary)*/
    /*No setting geom data as this is imported over from file, or is unnecessary for non-leaf cells*/

#if 0
    if((cell_level==0) &&(cell_num[0]=='0')) {
      printf("We're lookin 4 cell 0 into new cell with num %s\n",C->cell_num);
    }
#endif
    
    if(make_children(C) == ERROR)
      return(ERROR);

    establish_neighbour_after_refine(C -> children, ALL);

    if(assign_vtx_after_refine(C -> children, ALL, &Vtxlist, &Vtxlist_tail) == ERROR) /*Assign vertex pointers*/
      return(ERROR);
        
    for(i = 0; i < MAX_NUM_CHILDREN; i++) {
      if(add_to_leaf_list_front(C -> children[i], &(C -> Leaf_list_loc), &Leaves, &Leaves_tail) == ERROR) 
	return(ERROR);	      	      
    }

    if(C -> Leaf_list_loc != NULL) /*Only now delete the list node*/
      delete_from_leaf_list(C, &Leaves, &Leaves_tail);

    *this_cell = C;

#if 0
    if((cell_level==0) &&(cell_num[0]=='0')) {
      printf("We're just looked 4 cell 0 into new cell with num %s\n",C->cell_num);
    }
#endif
  }  

  return(NOERROR);
}
