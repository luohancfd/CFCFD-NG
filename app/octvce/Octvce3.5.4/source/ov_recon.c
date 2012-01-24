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

/**\file Source file for functions related to reconstruction, including gradient computation*/

#include <stdio.h>
#include "ov_kernel.h"
#include "ov_recon.h"
#include "ov_thermo.h"

extern Gas_dat gas_data;
extern short int twoD;
extern char many_limit;
extern int RK;
extern double switch_order;
extern double current_flow_time;

/*------------------------------------------------------------------*/

/**\brief Compute spatial gradients for all cells, particularly for merged cells where if normal cell can't be
   used as reconstruction matrix determinant 0, move onto small cells.  Otherwise small cells take same gradients as normal
   cell. If matrix determinant unavoidably zero, don't reconstruct here (set gradients to 0)*/

void get_grads_all_cells(List_leaf L, List_leaf Tail, short int new_cells_only)
{
  List_leaf node = L;
  List_merge merge_node;
  Cart_cell C, S;
  short int i, j;
  double det, grad[3][5];
  short int do_this_cell;
  
  while(node != NULL) {
    
    C = node -> cell_loc;

    do_this_cell = FALSE;
    /*Criterion for computing for these cells are - 
     (a) new_cells_only is FALSE and cell isn't small
     (b) new_cells_only is TRUE and cell is a new flow cell or cell is beside a new flow cell*/

#if DETON
    if(C -> un_det != TRUE)
#endif
      {
	if((new_cells_only == FALSE) && (C -> is_small_cell == FALSE))
	  do_this_cell = TRUE;
	else if(new_cells_only == TRUE) {
	  if(C -> just_created == TRUE)
	    do_this_cell = TRUE;
	  else if(C -> is_small_cell == FALSE) {
	    /*Now check and see if this cell is beside a newly created cell*/
	    for(i = 0; i < NUM_FACES; i++) {
	      if(C -> face_neighbours[i][0] != NULL) {
		if((C -> face_neighbours[i][1] == NULL) && (C -> face_neighbours[i][0] -> just_created == TRUE) && 
		   ((C -> flux_area[i][0] > 0) || (C -> flux_area[i][1] > 0) || (C->flux_area[i][2] > 0) || (C->flux_area[i][3] > 0))) {
		  do_this_cell = TRUE;
		  break;
		}
		else if(C -> face_neighbours[i][1] != NULL) {
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) {
		    if((C -> face_neighbours[i][j] -> just_created == TRUE) && (C -> flux_area[i][j] > 0)) {
		      do_this_cell = TRUE;
		      break;
		    }
		  }
		  if(do_this_cell == TRUE)
		    break;
		}
	      }
	    }
	  }
	}
      }

    if(do_this_cell == TRUE) {
            
      det = compute_spatial_grads(C, grad);

      if((det == 0) && (C -> Merge_data != NULL)) { /*Couldn't use any acceptable cells to get the spatial gradient - try small cells?*/
	merge_node = C -> Merge_data -> Linked_cells;

	while(merge_node != NULL) { 
	  S = merge_node -> cell_loc;

	  if(S -> is_small_cell == TRUE) {
	    det = compute_spatial_grads(S, grad);
	    if(det != 0) /*We could use this small cell to compute the gradients*/
	      break;
	  }
	  
	  merge_node = merge_node -> next;
	}	
      }

      if(det == 0) { /*No choice - no reconstruction will be employed for this cell/merged cluster*/
	for(i = 0; i < 3; i++)
	  for(j = 0; j < 5; j++)
	    grad[i][j] = 0;
      }
      
      if(twoD != 0) { /*As added precation, ensure all z gradient components are 0*/
	grad[2][0] = 0;
	grad[2][1] = 0;
	grad[2][2] = 0;
	grad[2][3] = 0;
	
	for(i = 0; i < 3; i++)
	  grad[i][4] = 0;	
      }

      /*Now finally let cells store this gradient*/
      if(C -> Merge_data != NULL) {
	merge_node = C -> Merge_data -> Linked_cells;
	
	while(merge_node != NULL) { /*All cells in merged cluster take same gradient values*/
	  S = merge_node -> cell_loc;

	  for(i = 0; i < 3; i++) {
	    S -> Flow_data.Grads.grad_rhof[i] = grad[i][0];
	    S -> Flow_data.Grads.grad_E[i] = grad[i][1];
	    S -> Flow_data.Grads.grad_u[i] = grad[i][2];
	    S -> Flow_data.Grads.grad_v[i] = grad[i][3];
	    S -> Flow_data.Grads.grad_w[i] = grad[i][4];
	  } 
	  
	  merge_node = merge_node -> next;
	}
      }
      else {
	for(i = 0; i < 3; i++) { 
	  C -> Flow_data.Grads.grad_rhof[i] = grad[i][0];
	  C -> Flow_data.Grads.grad_E[i] = grad[i][1];
	  C -> Flow_data.Grads.grad_u[i] = grad[i][2];
	  C -> Flow_data.Grads.grad_v[i] = grad[i][3];
	  C -> Flow_data.Grads.grad_w[i] = grad[i][4];	  
	}
      }
    }

    if(node == Tail)
      break;
    else node = node -> next;
  }
}

/*------------------------------------------------------------------*/

/**\brief Compute spatial gradients for a given cell using least-squares reconstruction, returns determinant
   of symmetric inverse reconstruction matrix*/

double compute_spatial_grads(Cart_cell C, double grad[][5])
{
  short int i, j, k, m;
  double dx, dy, dz, dxdy, dxdz, dydz, det;
  double recon_mat[3][3], R[3][5], ds[3], dW[5];

  dx = 0; dy = 0; dz = 0;
  dxdy = 0; dxdz = 0; dydz = 0;

  for(k = 0; k < 3; k++)  
    for(m = 0; m < 5; m++)  
      R[k][m] = 0;

  for(i = 0; i < NUM_FACES; i++) {
    if((C -> face_neighbours[i][0] != NULL) && (C -> face_neighbours[i][1] == NULL) &&
#if DETON
       (C -> face_neighbours[i][0] ->  un_det != TRUE) &&
#endif
       ((C->flux_area[i][0] > 0) || (C->flux_area[i][1] > 0) || (C->flux_area[i][2] > 0) || (C->flux_area[i][3] > 0)) &&
       ((C -> Merge_data == NULL) || (C -> Merge_data != C -> face_neighbours[i][0] -> Merge_data))) { /*1 cell in this direction*/ 

      for(k = 0; k < 3; k++)
	ds[k] = (C -> face_neighbours[i][0] -> centroid[k]) - (C -> centroid[k]);

      dx += SQR(ds[0]);
      dy += SQR(ds[1]);
      dz += SQR(ds[2]);
      dxdy += ds[0]*ds[1];
      dxdz += ds[0]*ds[2];
      dydz += ds[1]*ds[2];
	  
      dW[0] = C->face_neighbours[i][0]->Flow_data.State_vec.Rhof - C->Flow_data.State_vec.Rhof;
      dW[1] = (C->face_neighbours[i][0]->Flow_data.State_vec.RhoE)/(C->face_neighbours[i][0]->Flow_data.State_vec.Rhof) - 
	(C->Flow_data.State_vec.RhoE)/(C->Flow_data.State_vec.Rhof);
      dW[2] = (C->face_neighbours[i][0]->Flow_data.State_vec.RhoU)/(C->face_neighbours[i][0]->Flow_data.State_vec.Rhof) - 
	(C->Flow_data.State_vec.RhoU)/(C->Flow_data.State_vec.Rhof);
      dW[3] = (C->face_neighbours[i][0]->Flow_data.State_vec.RhoV)/(C->face_neighbours[i][0]->Flow_data.State_vec.Rhof) - 
	(C->Flow_data.State_vec.RhoV)/(C->Flow_data.State_vec.Rhof);
      dW[4] = (C->face_neighbours[i][0]->Flow_data.State_vec.RhoW)/(C->face_neighbours[i][0]->Flow_data.State_vec.Rhof) - 
	(C->Flow_data.State_vec.RhoW)/(C->Flow_data.State_vec.Rhof);
      
      for(k = 0; k < 3; k++)  /*Cycle through each direction*/
	for(m = 0; m < 5; m++)  /*Cycle through each flow variable*/
	  R[k][m] += ds[k]*dW[m];    	  
    }
    else if(C -> face_neighbours[i][1] != NULL) {
      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) {
	if((C -> flux_area[i][j] > 0) && 
#if DETON
	   (C -> face_neighbours[i][j] -> un_det != TRUE) &&
#endif
	   ((C -> Merge_data == NULL) || (C -> Merge_data != C -> face_neighbours[i][j] -> Merge_data))) { /*4 cells here*/

	  for(k = 0; k < 3; k++)
	    ds[k] = (C -> face_neighbours[i][j] -> centroid[k]) - (C -> centroid[k]);

	  dx += SQR(ds[0]);
	  dy += SQR(ds[1]);
	  dz += SQR(ds[2]);
	  dxdy += ds[0]*ds[1];
	  dxdz += ds[0]*ds[2];
	  dydz += ds[1]*ds[2];

	  dW[0] = C->face_neighbours[i][j]->Flow_data.State_vec.Rhof - C->Flow_data.State_vec.Rhof;
	  dW[1] = (C->face_neighbours[i][j]->Flow_data.State_vec.RhoE)/(C->face_neighbours[i][j]->Flow_data.State_vec.Rhof) - 
	    (C->Flow_data.State_vec.RhoE)/(C->Flow_data.State_vec.Rhof);
	  dW[2] = (C->face_neighbours[i][j]->Flow_data.State_vec.RhoU)/(C->face_neighbours[i][j]->Flow_data.State_vec.Rhof) - 
	    (C->Flow_data.State_vec.RhoU)/(C->Flow_data.State_vec.Rhof);
	  dW[3] = (C->face_neighbours[i][j]->Flow_data.State_vec.RhoV)/(C->face_neighbours[i][j]->Flow_data.State_vec.Rhof) - 
	    (C->Flow_data.State_vec.RhoV)/(C->Flow_data.State_vec.Rhof);
	  dW[4] = (C->face_neighbours[i][j]->Flow_data.State_vec.RhoW)/(C->face_neighbours[i][j]->Flow_data.State_vec.Rhof) - 
	    (C->Flow_data.State_vec.RhoW)/(C->Flow_data.State_vec.Rhof);
	  	  
	  for(k = 0; k < 3; k++)  
	    for(m = 0; m < 5; m++)  
	      R[k][m] += ds[k]*dW[m];
	}
      }
    }
  }

  if(twoD != 0) { /*We want to work with the 2x2 reconstruction matrix only, else 3x3 becomes singular*/
    det = dx*dy - SQR(dxdy);

    if(det == 0)
      return(det);

    recon_mat[0][0] = dy; recon_mat[0][1] = -dxdy;
    recon_mat[1][0] = -dxdy; recon_mat[1][1] = dx;

    for(i = 0; i < 2; i++)
      for(j = 0; j < 2; j++)
	recon_mat[i][j] = (recon_mat[i][j])/det;

    for(i = 0; i < 2; i++) {
      for(j = 0; j < 4; j++) {
	grad[i][j] = 0;
	for(k = 0; k < 2; k++)
	  grad[i][j] += recon_mat[i][k]*R[k][j];
      }
    }
  }
  else {
    /*Now we want to get the inverse of the inverse reconstruction matrix*/
    det = dx*(dy*dz - dydz*dydz) - dxdy*(dxdy*dz - dydz*dxdz) + dxdz*(dxdy*dydz - dxdz*dy);

    if(det == 0)
      return(det);

    /*Determinant would never be 0 or else cell won't be on the list*/
    recon_mat[0][0] = dy*dz - dydz*dydz; recon_mat[0][1] = dxdz*dydz - dxdy*dz; recon_mat[0][2] = dxdy*dydz - dy*dxdz;
    recon_mat[1][0] = dxdz*dydz - dxdy*dz; recon_mat[1][1] = dx*dz - dxdz*dxdz; recon_mat[1][2] = dxdy*dxdz - dx*dydz;
    recon_mat[2][0] = dxdy*dydz - dxdz*dy; recon_mat[2][1] = dxdz*dxdy - dx*dydz; recon_mat[2][2] = dx*dy - dxdy*dxdy;

    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	recon_mat[i][j] = (recon_mat[i][j])/det; /*Matrix finally obtained*/

    /*Now multiply reconstruction matrix with matrix of 'right hand' quantities R*/

    for(i = 0; i < 3; i++) { /*For each row of the gradient matrix*/
      for(j = 0; j < 5; j++) { /*For each column of the gradient matrix*/
	grad[i][j] = 0;
	for(k = 0; k < 3; k++) /*Go through 'inner' dimension*/
	  grad[i][j] += recon_mat[i][k]*R[k][j];
      }
    }
  }
 
  return(det);
}

/*------------------------------------------------------------------*/

/**\brief Compute gradient limiter for each cell (of the min-mod variety).  Small cells will have the same limite value as
   their normal counterparts*/

void compute_limiter(List_leaf L, List_leaf Tail, short int new_cells_only)
{
  List_leaf node = L;
  List_merge merge_node;
  Cart_cell C, S;
  short int i, j, k, do_this_cell;
  double dist, phi, phimin, q, qcen, qmax, qmin, lims[5], minlim;
  double dr[3], grad[3];
  double num_var, num_vtx, num_faces;

  if(twoD != 0) {
    num_var = 4;
    num_vtx = 4;
    num_faces = 4;
  }
  else {
    num_var = 5;
    num_vtx = 8;
    num_faces = 6;
  }
  q = 0; qcen = 0; minlim = 0;
  while(node != NULL) {
    
    C = node -> cell_loc;

    do_this_cell = FALSE;
    /*Criterion for computing for these cells are - 
     (a) new_cells_only is FALSE and cell isn't small
     (b) new_cells_only is TRUE and cell is a new flow cell or cell is beside a new flow cell*/

#if DETON
    if(C -> un_det != TRUE)
#endif
      {
	if((new_cells_only == FALSE) && (C -> is_small_cell == FALSE))
	  do_this_cell = TRUE;
	else if(new_cells_only == TRUE) {
	  if(C -> just_created == TRUE)
	    do_this_cell = TRUE;
	  else if(C -> is_small_cell == FALSE) {
	    /*Now check and see if this cell is beside a newly created cell*/
	    for(i = 0; i < NUM_FACES; i++) {
	      if(C -> face_neighbours[i][0] != NULL) {
		if((C -> face_neighbours[i][1] == NULL) && (C -> face_neighbours[i][0] -> just_created == TRUE) && 
		   ((C -> flux_area[i][0] > 0) || (C -> flux_area[i][1] > 0) || (C->flux_area[i][2] > 0) || (C->flux_area[i][3] > 0))) {
		  do_this_cell = TRUE;
		  break;
		}
		else if(C -> face_neighbours[i][1] != NULL) {
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) {
		    if((C -> face_neighbours[i][j] -> just_created == TRUE) && (C -> flux_area[i][j] > 0)) {
		      do_this_cell = TRUE;
		      break;
		    }
		  }
		  if(do_this_cell == TRUE)
		    break;
		}
	      }
	    }
	  }
	}
      }

    if(do_this_cell == TRUE) {

      dist = 0.5*CALC_CELL_EDGE_LENGTH(C); /*Half a cell length*/

      C -> Flow_data.Lim.rhof_lim = 0;
      C -> Flow_data.Lim.E_lim = 0;
      C -> Flow_data.Lim.u_lim = 0;
      C -> Flow_data.Lim.v_lim = 0;
      C -> Flow_data.Lim.w_lim = 0;
      C -> Flow_data.Lim.lim = 0;

      for(i = 0; i < num_var; i++) { /*Go through all flow variables (limiter values for each)*/

	/*If only 1 species, 5th limiter value computed redundantly*/

	switch(i) {
	case 0:
	  qcen = C -> Flow_data.State_vec.Rhof; break;
	case 1:
	  qcen = (C -> Flow_data.State_vec.RhoE)/(C -> Flow_data.State_vec.Rhof); break;
	case 2:
	  qcen = (C -> Flow_data.State_vec.RhoU)/(C -> Flow_data.State_vec.Rhof); break;
	case 3:
	  qcen = (C -> Flow_data.State_vec.RhoV)/(C -> Flow_data.State_vec.Rhof); break;
	case 4:
	  qcen = (C -> Flow_data.State_vec.RhoW)/(C -> Flow_data.State_vec.Rhof); break;	
	}

	qmin = qcen; 
	qmax = qcen;

	/*Get minimum and maximum quantities*/

	for(j = 0; j < num_faces; j++) {
	  if((C -> face_neighbours[j][0] != NULL) && (C -> face_neighbours[j][1] == NULL) &&
#if DETON
	     (C -> face_neighbours[j][0] -> un_det != TRUE) &&
#endif
	     ((C->flux_area[j][0] > 0) || (C->flux_area[j][1] > 0) || (C->flux_area[j][2] > 0) || (C->flux_area[j][3] > 0)) &&
	     ((C -> Merge_data == NULL) || (C -> Merge_data != C -> face_neighbours[j][0] -> Merge_data))) {

	    switch(i) {
	    case 0:
	      q = C -> face_neighbours[j][0] -> Flow_data.State_vec.Rhof; break;
	    case 1:
	      q = (C -> face_neighbours[j][0] -> Flow_data.State_vec.RhoE)/(C -> face_neighbours[j][0] -> Flow_data.State_vec.Rhof); break;
	    case 2:
	      q = (C -> face_neighbours[j][0] -> Flow_data.State_vec.RhoU)/(C -> face_neighbours[j][0] -> Flow_data.State_vec.Rhof); break;
	    case 3:
	      q = (C -> face_neighbours[j][0] -> Flow_data.State_vec.RhoV)/(C -> face_neighbours[j][0] -> Flow_data.State_vec.Rhof); break;
	    case 4:
	      q = (C -> face_neighbours[j][0] -> Flow_data.State_vec.RhoW)/(C -> face_neighbours[j][0] -> Flow_data.State_vec.Rhof); break;
	    }

	    if(q < qmin)
	      qmin = q;
	    if(q > qmax)
	      qmax = q;
	  }
	  else if(C -> face_neighbours[j][1] != NULL) {
	    for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) {
	      if((C -> flux_area[j][k] > 0) && 
#if DETON
		 (C -> face_neighbours[j][k] -> un_det != TRUE) &&
#endif
		 ((C -> Merge_data == NULL) || (C -> Merge_data != C -> face_neighbours[j][k] -> Merge_data))) { /*4 cells here*/
		switch(i) {
		case 0:
		  q = C -> face_neighbours[j][k] -> Flow_data.State_vec.Rhof; break;
		case 1:
		  q = (C->face_neighbours[j][k]->Flow_data.State_vec.RhoE)/(C->face_neighbours[j][k]->Flow_data.State_vec.Rhof); break;
		case 2:
		  q = (C->face_neighbours[j][k]->Flow_data.State_vec.RhoU)/(C->face_neighbours[j][k]->Flow_data.State_vec.Rhof); break;
		case 3:
		  q = (C->face_neighbours[j][k]->Flow_data.State_vec.RhoV)/(C->face_neighbours[j][k]->Flow_data.State_vec.Rhof); break;
		case 4:
		  q = (C->face_neighbours[j][k]->Flow_data.State_vec.RhoW)/(C->face_neighbours[j][k]->Flow_data.State_vec.Rhof); break;
		}

		if(q < qmin)
		  qmin = q;
		if(q > qmax)
		  qmax = q;
	      }
	    }
	  }
	} 

	switch(i) {
	case 0:
	  for(k = 0; k < 3; k++)
	    grad[k] = C -> Flow_data.Grads.grad_rhof[k];
	  break;
	case 1:
	  for(k = 0; k < 3; k++)
	    grad[k] = C -> Flow_data.Grads.grad_E[k];
	  break;
	case 2:
	  for(k = 0; k < 3; k++)
	    grad[k] = C -> Flow_data.Grads.grad_u[k];
	  break;
	case 3:
	  for(k = 0; k < 3; k++)
	    grad[k] = C -> Flow_data.Grads.grad_v[k];
	  break;
	case 4:
	  for(k = 0; k < 3; k++)
	    grad[k] = C -> Flow_data.Grads.grad_w[k];
	  break;
	}

	phimin = 1;

	for(j = 0; j < num_vtx; j++) { /*Go through all verticies*/
	  switch(j) {
	  case 0: /*Vertex 0*/
	    dr[0] = -1; dr[1] = -1; dr[2] = -1; break;
	  case 1:
	    dr[0] = -1; dr[1] = 1; dr[2] = -1; break;
	  case 2:
	    dr[0] = 1; dr[1] = -1; dr[2] = -1; break;
	  case 3:
	    dr[0] = 1; dr[1] = 1; dr[2] = -1; break;
	  case 4:
	    dr[0] = -1; dr[1] = -1; dr[2] = 1; break;
	  case 5:
	    dr[0] = -1; dr[1] = 1; dr[2] = 1; break;
	  case 6:
	    dr[0] = 1; dr[1] = -1; dr[2] = 1; break;
	  case 7:
	    dr[0] = 1; dr[1] = 1; dr[2] = 1; break;
	  }

	  for(k = 0; k < 3; k++)
	    dr[k] = dist*dr[k]; /*Distance to vertex calculated*/ 

	  if(twoD != 0)
	    q = qcen + grad[0]*dr[0] + grad[1]*dr[1];
	  else q = qcen + grad[0]*dr[0] + grad[1]*dr[1] + grad[2]*dr[2]; /*Extrapolated value without limiting calculated*/ 
	
	  /*Compute phi for each vertex*/
	  if(q - qcen > 0)
	    phi = MIN(1, (qmax - qcen)/(q - qcen));
	  else if(q - qcen < 0)
	    phi = MIN(1, (qmin - qcen)/(q - qcen));
	  else phi = 1;

	  if(phi < phimin)
	    phimin = phi; /*Get minimum phi over all verticies (which will be the limiter values)*/
	}

	if(many_limit == 'n')
	  lims[i] = phimin; /*Get limiter for each value before finding minimum out of them all*/
        else {
	  switch(i) {
	  case 0:
	    C -> Flow_data.Lim.rhof_lim = phimin; break;
	  case 1:
	    C -> Flow_data.Lim.E_lim = phimin; break;
	  case 2:
	    C -> Flow_data.Lim.u_lim = phimin; break;
	  case 3:
	    C -> Flow_data.Lim.v_lim = phimin; break;
	  case 4:
	    C -> Flow_data.Lim.w_lim = phimin; break;
	  }
	}

#if NOTINCLUDE
	if(phimin > 0)
	  printf("PHIMIN > 0, it's %g\n",phimin);
#endif

      } /*Finished going through all flow variables*/


      if(many_limit == 'n') {
	/*Get total minimum limiter*/
	minlim = lims[0];
	for(i = 1; i < num_var; i++)
	  if(lims[i] < minlim)
	    minlim = lims[i];
      
	C -> Flow_data.Lim.lim = minlim;
      }

      if(C -> Merge_data != NULL) {
	merge_node = C -> Merge_data -> Linked_cells;

	while(merge_node != NULL) { /*Merged cells take same limiter values as their normal counterparts*/
	  S = merge_node -> cell_loc;

	  if(many_limit == 'n') {
	    if(S -> is_small_cell == TRUE) 
	      S -> Flow_data.Lim.lim = minlim;
	  }
	  else {
	    if(S -> is_small_cell == TRUE) {
	      S -> Flow_data.Lim.rhof_lim = C -> Flow_data.Lim.rhof_lim;
	      S -> Flow_data.Lim.E_lim = C -> Flow_data.Lim.E_lim;
	      S -> Flow_data.Lim.u_lim = C -> Flow_data.Lim.u_lim;
	      S -> Flow_data.Lim.v_lim = C -> Flow_data.Lim.v_lim;
	      S -> Flow_data.Lim.w_lim = C -> Flow_data.Lim.w_lim;
	    }
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

/*------------------------------------------------------------------*/

/**\brief Reconstruct state vectors for children when parent is refined.  This must be done conservatively - thus
   all children have the same cell-centered state as their parent*/

void reconstruct_refined_children(Cart_cell child_cells[])
{
  short int i;
  Cart_cell Parent = child_cells[0] -> parent;
  
#if 0
  for(i = 0; i < MAX_NUM_CHILDREN; i++)
    {
      if(child_cells[i] -> cell_type != SOLID)
	{
	  child_cells[i] -> Flow_data.State_vec.Rhof = Parent -> Flow_data.State_vec.Rhof;
	  child_cells[i] -> Flow_data.State_vec.RhoU = Parent -> Flow_data.State_vec.RhoU;
	  child_cells[i] -> Flow_data.State_vec.RhoV = Parent -> Flow_data.State_vec.RhoV;
	  child_cells[i] -> Flow_data.State_vec.RhoW = Parent -> Flow_data.State_vec.RhoW;
	  child_cells[i] -> Flow_data.State_vec.RhoE = Parent -> Flow_data.State_vec.RhoE;
	  child_cells[i] -> Flow_data.State_vec.Pres = Parent -> Flow_data.State_vec.Pres;
	  child_cells[i] -> Flow_data.State_vec.T = Parent -> Flow_data.State_vec.T;
	  
	  if(gas_data.num_species > 1) {
	    child_cells[i] -> Flow_data.State_vec.rhof_prod = Parent -> Flow_data.State_vec.rhof_prod;
#if AFTERBURN
	    child_cells[i] -> Flow_data.State_vec.rhof_prod_ub = Parent -> Flow_data.State_vec.rhof_prod_ub;
#endif
	  }
	}
    }
#else
  /*Do reconstruction, but ensure conservativity by first reconstructing then subtracting the 'error'*/

  short int num_children;
  double Rhofsum, RhoUsum, RhoVsum, RhoWsum, RhoEsum, Rhof[8], RhoU[8], RhoV[8], RhoW[8], RhoE[8], prod_m_frac;
  double rhof_prod[8], vol;
  double dRhof, dRhoU, dRhoV, dRhoW, dRhoE, R;

  R = (gas_data.Cv_amb)*((gas_data.gam_amb) - 1.0);

#if AFTERBURN
  double rhof_prod_ub[8];
  double prod_m_frac_ub;
#endif
  prod_m_frac = 0;
  if(twoD != 0)
    num_children = 4;
  else num_children = 8;
  
  Rhofsum = 0; RhoUsum = 0; RhoVsum = 0; RhoWsum = 0; RhoEsum = 0, vol = 0; 

  if(gas_data.num_species > 1) {
    prod_m_frac = (Parent -> Flow_data.State_vec.rhof_prod)/(Parent -> Flow_data.State_vec.Rhof);
#if AFTERBURN
    prod_m_frac_ub = (Parent -> Flow_data.State_vec.rhof_prod_ub)/(Parent -> Flow_data.State_vec.rhof_prod); 
#endif
  }

  for(i = 0; i < num_children; i++) {
    if(child_cells[i] -> cell_type != SOLID) {
      reconstruct_flow(Parent, IRRELEVANT, IRRELEVANT, i, 0.5, 'R', &(Rhof[i]), NULL, 'n');
      reconstruct_flow(Parent, IRRELEVANT, IRRELEVANT, i, 0.5, 'U', &(RhoU[i]), NULL, 'n');
      reconstruct_flow(Parent, IRRELEVANT, IRRELEVANT, i, 0.5, 'V', &(RhoV[i]), NULL, 'n');
      reconstruct_flow(Parent, IRRELEVANT, IRRELEVANT, i, 0.5, 'W', &(RhoW[i]), NULL, 'n');
      reconstruct_flow(Parent, IRRELEVANT, IRRELEVANT, i, 0.5, 'E', &(RhoE[i]), NULL, 'n');

      RhoU[i] = RhoU[i]*Rhof[i]; /*It was just plain velocity*/
      RhoV[i] = RhoV[i]*Rhof[i];
      RhoW[i] = RhoW[i]*Rhof[i];
      RhoE[i] = RhoE[i]*Rhof[i];

      Rhofsum += Rhof[i]*(child_cells[i]->cell_volume); /*Get total conserved quantity*/
      RhoUsum += RhoU[i]*(child_cells[i]->cell_volume);
      RhoVsum += RhoV[i]*(child_cells[i]->cell_volume);
      RhoWsum += RhoW[i]*(child_cells[i]->cell_volume);
      RhoEsum += RhoE[i]*(child_cells[i]->cell_volume);

      vol += (child_cells[i]->cell_volume);
    }
  }

  /*Get difference per cell (to subtract from each cell to ensure conservativity)*/
  dRhof = (vol*(Parent -> Flow_data.State_vec.Rhof) - Rhofsum)/((double) i);
  dRhoU = (vol*(Parent -> Flow_data.State_vec.RhoU) - RhoUsum)/((double) i);
  dRhoV = (vol*(Parent -> Flow_data.State_vec.RhoV) - RhoVsum)/((double) i);
  dRhoW = (vol*(Parent -> Flow_data.State_vec.RhoW) - RhoWsum)/((double) i);
  dRhoE = (vol*(Parent -> Flow_data.State_vec.RhoE) - RhoEsum)/((double) i);
  
  for(i = 0; i < num_children; i++) {
    if(child_cells[i] -> cell_type != SOLID) {
      Rhof[i] = Rhof[i] + dRhof/(child_cells[i]->cell_volume); /*Finally add/subtract to get conserved value*/
      RhoU[i] = RhoU[i] + dRhoU/(child_cells[i]->cell_volume);
      RhoV[i] = RhoV[i] + dRhoV/(child_cells[i]->cell_volume);
      RhoW[i] = RhoW[i] + dRhoW/(child_cells[i]->cell_volume);
      RhoE[i] = RhoE[i] + dRhoE/(child_cells[i]->cell_volume);

      child_cells[i] -> Flow_data.State_vec.Rhof = Rhof[i];
      child_cells[i] -> Flow_data.State_vec.RhoU = RhoU[i];
      child_cells[i] -> Flow_data.State_vec.RhoV = RhoV[i];
      child_cells[i] -> Flow_data.State_vec.RhoW = RhoW[i];
      child_cells[i] -> Flow_data.State_vec.RhoE = RhoE[i];

      if(gas_data.num_species > 1) {
	rhof_prod[i] = prod_m_frac*Rhof[i]; /*Mass fractions constant though density can vary*/
	child_cells[i] -> Flow_data.State_vec.rhof_prod = rhof_prod[i];
#if AFTERBURN
	rhof_prod_ub[i] = prod_m_frac_ub*rhof_prod[i]; 
	child_cells[i] -> Flow_data.State_vec.rhof_prod_ub = rhof_prod_ub[i];
#endif

	get_mixture_Pres(&(child_cells[i] -> Flow_data.State_vec));
      }
      else {
	child_cells[i] -> Flow_data.State_vec.Pres = ((gas_data.gam_amb)-1.0)*
	  ((child_cells[i]->Flow_data.State_vec.RhoE) - 
	   0.5*(SQR(child_cells[i]->Flow_data.State_vec.RhoU) + SQR(child_cells[i]->Flow_data.State_vec.RhoV) +
		SQR(child_cells[i]->Flow_data.State_vec.RhoW))/(child_cells[i]->Flow_data.State_vec.Rhof));
	child_cells[i] -> Flow_data.State_vec.T = (child_cells[i] -> Flow_data.State_vec.Pres)/
	  (R*(child_cells[i]->Flow_data.State_vec.Rhof));
      }      
    }
  }
#endif
}

/*------------------------------------------------------------------*/

/**\brief Reconstruct state vectors for cell when it is coarsened based on children.  To ensure conservativity, take
   volume-weighted average of childrens' state vectors (for conserved quantities, of course).*/

void reconstruct_coarsened_cell(Cart_cell C)
{
  double Rhofsum, RhoUsum, RhoVsum, RhoWsum, RhoEsum, rhof_prod_sum, child_vol, vol;
  short int i;

#if AFTERBURN
  double rhof_prod_ub_sum = 0;
#endif
  
  Rhofsum = 0; RhoUsum = 0; RhoVsum = 0; RhoWsum = 0; RhoEsum = 0; rhof_prod_sum = 0;
  vol = C -> cell_volume;

  for(i = 0; i < MAX_NUM_CHILDREN; i++)
    {
      if(C -> children[i] -> cell_type != SOLID) /*Only use computational cells*/
	{
	  child_vol = C -> children[i] -> cell_volume;

	  Rhofsum += (C -> children[i] -> Flow_data.State_vec.Rhof)*child_vol;
	  RhoUsum += (C -> children[i] -> Flow_data.State_vec.RhoU)*child_vol;
	  RhoVsum += (C -> children[i] -> Flow_data.State_vec.RhoV)*child_vol;
	  RhoWsum += (C -> children[i] -> Flow_data.State_vec.RhoW)*child_vol;
	  RhoEsum += (C -> children[i] -> Flow_data.State_vec.RhoE)*child_vol;

	  if(gas_data.num_species > 1) {
	    rhof_prod_sum += (C -> children[i] -> Flow_data.State_vec.rhof_prod)*child_vol;
#if AFTERBURN
	    rhof_prod_ub_sum += (C -> children[i] -> Flow_data.State_vec.rhof_prod_ub)*child_vol;
#endif
	  }
	}
    }
  
  C -> Flow_data.State_vec.Rhof = Rhofsum/vol;
  C -> Flow_data.State_vec.RhoU = RhoUsum/vol;
  C -> Flow_data.State_vec.RhoV = RhoVsum/vol;
  C -> Flow_data.State_vec.RhoW = RhoWsum/vol;
  C -> Flow_data.State_vec.RhoE = RhoEsum/vol;
    
  if(gas_data.num_species > 1) {
    C -> Flow_data.State_vec.rhof_prod = rhof_prod_sum/vol;
#if AFTERBURN
    C -> Flow_data.State_vec.rhof_prod = rhof_prod_ub_sum/vol;
#endif

    get_mixture_Pres(&(C -> Flow_data.State_vec)); /*For EOS consistency - will get temperature and pressure*/
  }
  else C -> Flow_data.State_vec.Pres = ((gas_data.gam_amb)-1.0)*
	 ((C->Flow_data.State_vec.RhoE) - 0.5*(SQR(C->Flow_data.State_vec.RhoU) + SQR(C->Flow_data.State_vec.RhoV) +
					       SQR(C->Flow_data.State_vec.RhoW))/(C->Flow_data.State_vec.Rhof));    
}

/*------------------------------------------------------------------*/

/**\brief Reconstruct flow variables at each cell interface using min-mod style limiting*/

State_vector reconstruct_flow(Cart_cell C, short int dir, short int quad, short int vtx, short int distance, int var, 
			      double *q, double *loc, char full_xtrp)
{
  State_vector sv;
  short int i;
  double dr[3];  
  double dist = 0.5*CALC_CELL_EDGE_LENGTH(C); /*Half a cell length*/
  double lim, rhof_lim, E_lim, u_lim, v_lim, w_lim;
  lim = 0; rhof_lim = 0; E_lim = 0; u_lim = 0; v_lim = 0; w_lim = 0;
  sv.rhof_prod = 0;

  if((switch_order > 0) && (current_flow_time < switch_order)) {
    lim = 0.0;
    rhof_lim = 0.0; E_lim = 0.0; u_lim = 0.0; v_lim = 0.0; w_lim = 0.0;
  }
  else {
    if(full_xtrp == 'y') {
      lim = 1.0; 
      rhof_lim = 1.0; E_lim = 1.0; u_lim = 1.0; v_lim = 1.0; w_lim = 1.0;
    }
    else {
      if(
#if 1
C -> cell_type == INTERSECTED
#else 
0
#endif
) { /*Don't reconstruct for intersected cells*/
	lim = 0.0;
	rhof_lim = 0.0; E_lim = 0.0; u_lim = 0.0; v_lim = 0.0; w_lim = 0.0;
      }
      else {
	if(many_limit == 'n')
	  lim = C -> Flow_data.Lim.lim;
	else {
	  rhof_lim = C -> Flow_data.Lim.rhof_lim;
	  E_lim = C -> Flow_data.Lim.E_lim;
	  u_lim = C -> Flow_data.Lim.u_lim;
	  v_lim = C -> Flow_data.Lim.v_lim;
	  w_lim = C -> Flow_data.Lim.w_lim;
	}
      }
    }
  }

#if !NOTINCLUDE
       
  if(dir != IRRELEVANT) { /*Reconstruct to interface*/
    switch(dir) {
    case EST: 
      dr[0] = 1; dr[1] = 0; dr[2] = 0; break;
    case WST:
      dr[0] = -1; dr[1] = 0; dr[2] = 0; break;
    case NTH:
      dr[0] = 0; dr[1] = 1; dr[2] = 0; break;
    case STH:
      dr[0] = 0; dr[1] = -1; dr[2] = 0; break;
    case UPR:
      dr[0] = 0; dr[1] = 0; dr[2] = 1; break;
    case LWR:
      dr[0] = 0; dr[1] = 0; dr[2] = -1; break;  
    }

    if(quad != IRRELEVANT) { /*4 neighbours on this face*/
      if((dir == UPR) || (dir == LWR)) {
	switch(quad) {
	case 0:
	  dr[0] += -0.5; dr[1] += -0.5; break; /*Need to go 0.5 of 0.5 cell length to interface centroid*/
	case 1:
	  dr[0] += 0.5; dr[1] += 0.5; break;
	case 2:
	  dr[0] += 0.5; dr[1] += -0.5; break;
	case 3:
	  dr[0] += 0.5; dr[1] += 0.5; break;
	}
      }
      else if((dir == NTH) || (dir == STH)) {
	switch(quad) {
	case 0:
	  dr[0] += -0.5; dr[2] += -0.5; break;
	case 1:
	  dr[0] += 0.5; dr[2] += -0.5; break;
	case 2:
	  dr[0] += -0.5; dr[2] += 0.5; break;
	case 3:
	  dr[0] += 0.5; dr[2] += 0.5; break;
	}
      }
      else if((dir == EST) || (dir == WST)) {
	switch(quad) {
	case 0:
	  dr[1] += -0.5; dr[2] += -0.5; break;
	case 1:
	  dr[1] += 0.5; dr[2] += -0.5; break;
	case 2:
	  dr[1] += -0.5; dr[2] += 0.5; break;
	case 3:
	  dr[1] += 0.5; dr[2] += 0.5; break;
	}
      }
    }

    for(i = 0; i < 3; i++) /*Scale to cell size*/
      dr[i] = dist*dr[i];
  }
  else if(vtx != IRRELEVANT) { /*Reconstruct to cell corner*/
    switch(vtx) {
    case 0: /*Vertex 0*/
      dr[0] = -1; dr[1] = -1; dr[2] = -1; break;
    case 1:
      dr[0] = -1; dr[1] = 1; dr[2] = -1; break;
    case 2:
      dr[0] = 1; dr[1] = -1; dr[2] = -1; break;
    case 3:
      dr[0] = 1; dr[1] = 1; dr[2] = -1; break;
    case 4:
      dr[0] = -1; dr[1] = -1; dr[2] = 1; break;
    case 5:
      dr[0] = -1; dr[1] = 1; dr[2] = 1; break;
    case 6:
      dr[0] = 1; dr[1] = -1; dr[2] = 1; break;
    case 7:
      dr[0] = 1; dr[1] = 1; dr[2] = 1; break;
    }

    for(i = 0; i < 3; i++)
      dr[i] = dist*distance*dr[i]; /*May wish to only reconstruct half way to cell corner e.g. refining cell*/
  }
  else { /*Want to reconstruct to an arbitrary location*/
    dr[0] = loc[0] - C->centroid[0];
    dr[1] = loc[1] - C->centroid[1];
    dr[2] = loc[2] - C->centroid[2];
  }

  if(twoD != 0) 
    dr[2] = 0; /*Z-coordinate is irrelevant*/
  
  if(var == EOF) {
    if(many_limit == 'n') {
      sv.Rhof = C->Flow_data.State_vec.Rhof + lim*
	((C->Flow_data.Grads.grad_rhof[0])*dr[0] + (C->Flow_data.Grads.grad_rhof[1])*dr[1] + (C->Flow_data.Grads.grad_rhof[2])*dr[2]);
      sv.RhoE = (sv.Rhof)*
	((C->Flow_data.State_vec.RhoE)/(C->Flow_data.State_vec.Rhof) + lim*
	 ((C->Flow_data.Grads.grad_E[0])*dr[0] + (C->Flow_data.Grads.grad_E[1])*dr[1] + (C->Flow_data.Grads.grad_E[2])*dr[2]));
      sv.RhoU = (sv.Rhof)*
	((C->Flow_data.State_vec.RhoU)/(C->Flow_data.State_vec.Rhof) + lim*
	 ((C->Flow_data.Grads.grad_u[0])*dr[0] + (C->Flow_data.Grads.grad_u[1])*dr[1] + (C->Flow_data.Grads.grad_u[2])*dr[2]));
      sv.RhoV = (sv.Rhof)*
	((C->Flow_data.State_vec.RhoV)/(C->Flow_data.State_vec.Rhof) + lim*
	 ((C->Flow_data.Grads.grad_v[0])*dr[0] + (C->Flow_data.Grads.grad_v[1])*dr[1] + (C->Flow_data.Grads.grad_v[2])*dr[2]));
      sv.RhoW = (sv.Rhof)*
	((C->Flow_data.State_vec.RhoW)/(C->Flow_data.State_vec.Rhof) + lim*
	 ((C->Flow_data.Grads.grad_w[0])*dr[0] + (C->Flow_data.Grads.grad_w[1])*dr[1] + (C->Flow_data.Grads.grad_w[2])*dr[2]));
    }
    else {
      sv.Rhof = C->Flow_data.State_vec.Rhof + rhof_lim*
	((C->Flow_data.Grads.grad_rhof[0])*dr[0] + (C->Flow_data.Grads.grad_rhof[1])*dr[1] + (C->Flow_data.Grads.grad_rhof[2])*dr[2]);
      sv.RhoE = (sv.Rhof)*
	((C->Flow_data.State_vec.RhoE)/(C->Flow_data.State_vec.Rhof) + E_lim*
	 ((C->Flow_data.Grads.grad_E[0])*dr[0] + (C->Flow_data.Grads.grad_E[1])*dr[1] + (C->Flow_data.Grads.grad_E[2])*dr[2]));
      sv.RhoU = (sv.Rhof)*
	((C->Flow_data.State_vec.RhoU)/(C->Flow_data.State_vec.Rhof) + u_lim*
	 ((C->Flow_data.Grads.grad_u[0])*dr[0] + (C->Flow_data.Grads.grad_u[1])*dr[1] + (C->Flow_data.Grads.grad_u[2])*dr[2]));
      sv.RhoV = (sv.Rhof)*
	((C->Flow_data.State_vec.RhoV)/(C->Flow_data.State_vec.Rhof) + v_lim*
	 ((C->Flow_data.Grads.grad_v[0])*dr[0] + (C->Flow_data.Grads.grad_v[1])*dr[1] + (C->Flow_data.Grads.grad_v[2])*dr[2]));
      sv.RhoW = (sv.Rhof)*
	((C->Flow_data.State_vec.RhoW)/(C->Flow_data.State_vec.Rhof) + w_lim*
	 ((C->Flow_data.Grads.grad_w[0])*dr[0] + (C->Flow_data.Grads.grad_w[1])*dr[1] + (C->Flow_data.Grads.grad_w[2])*dr[2]));
    }

    if((sv.RhoE - 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof)) <= 0) {
      /*When gradients are very high could this happen? - then don't reconstruct*/

#if 1
      sv.Rhof = C -> Flow_data.State_vec.Rhof;
      sv.RhoU = C -> Flow_data.State_vec.RhoU;
      sv.RhoV = C -> Flow_data.State_vec.RhoV;
      sv.RhoW = C -> Flow_data.State_vec.RhoW;
      sv.RhoE = C -> Flow_data.State_vec.RhoE;
      sv.Pres = C -> Flow_data.State_vec.Pres;

      if(gas_data.num_species > 1) {
	sv.rhof_prod = C -> Flow_data.State_vec.rhof_prod;
#if AFTERBURN
	sv.rhof_prod_ub = C -> Flow_data.State_vec.rhof_prod_ub;
#endif
      }
#endif
    }
    else {
      if(gas_data.num_species > 1) {

	sv.rhof_prod = ((C -> Flow_data.State_vec.rhof_prod)/(C -> Flow_data.State_vec.Rhof))*(sv.Rhof);
	/*Mass fractions are constant throughout cell though density can vary*/

#if AFTERBURN      
	sv.rhof_prod_ub = ((C -> Flow_data.State_vec.rhof_prod_ub)/(C -> Flow_data.State_vec.rhof_prod))*(sv.rhof_prod);
#endif

	get_mixture_Pres(&sv); /*Let's hope this is a good interface value*/
      }
      else sv.Pres = ((sv.RhoE) - 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof))*((gas_data.gam_amb) - 1.0); 
    }
  }    
  else { /*Want to reconstruct a specific variable*/
    if(many_limit == 'n') {
      if(var == 'P') /*Only want to reconstruct pressure*/
	*q = C->Flow_data.State_vec.Pres; /*Don't bother using EOS (as it's expensive) - for I/O, cell-centered would do methinks*/
      else if(var == 'R')
	*q = C->Flow_data.State_vec.Rhof + (lim)*
	  ((C->Flow_data.Grads.grad_rhof[0])*dr[0] + (C->Flow_data.Grads.grad_rhof[1])*dr[1] + (C->Flow_data.Grads.grad_rhof[2])*dr[2]);
      else if(var == 'U')
	*q = (C->Flow_data.State_vec.RhoU)/(C->Flow_data.State_vec.Rhof) + (lim)*
	  ((C->Flow_data.Grads.grad_u[0])*dr[0] + (C->Flow_data.Grads.grad_u[1])*dr[1] + (C->Flow_data.Grads.grad_u[2])*dr[2]);
      else if(var == 'V')
	*q = (C->Flow_data.State_vec.RhoV)/(C->Flow_data.State_vec.Rhof) + (lim)*
	  ((C->Flow_data.Grads.grad_v[0])*dr[0] + (C->Flow_data.Grads.grad_v[1])*dr[1] + (C->Flow_data.Grads.grad_v[2])*dr[2]);
      else if(var == 'W')
	*q = (C->Flow_data.State_vec.RhoW)/(C->Flow_data.State_vec.Rhof) + (lim)*
	  ((C->Flow_data.Grads.grad_w[0])*dr[0] + (C->Flow_data.Grads.grad_w[1])*dr[1] + (C->Flow_data.Grads.grad_w[2])*dr[2]);
      else if(var == 'E')
	*q = (C->Flow_data.State_vec.RhoE)/(C->Flow_data.State_vec.Rhof) + (lim)*
	  ((C->Flow_data.Grads.grad_E[0])*dr[0] + (C->Flow_data.Grads.grad_E[1])*dr[1] + (C->Flow_data.Grads.grad_E[2])*dr[2]);
    }
    else {
      if(var == 'P') /*Only want to reconstruct pressure*/
	*q = C->Flow_data.State_vec.Pres;
      else if(var == 'R')
	*q = C->Flow_data.State_vec.Rhof + (rhof_lim)*
	  ((C->Flow_data.Grads.grad_rhof[0])*dr[0] + (C->Flow_data.Grads.grad_rhof[1])*dr[1] + (C->Flow_data.Grads.grad_rhof[2])*dr[2]);
      else if(var == 'U')
	*q = (C->Flow_data.State_vec.RhoU)/(C->Flow_data.State_vec.Rhof) + (u_lim)*
	  ((C->Flow_data.Grads.grad_u[0])*dr[0] + (C->Flow_data.Grads.grad_u[1])*dr[1] + (C->Flow_data.Grads.grad_u[2])*dr[2]);
      else if(var == 'V')
	*q = (C->Flow_data.State_vec.RhoV)/(C->Flow_data.State_vec.Rhof) + (v_lim)*
	  ((C->Flow_data.Grads.grad_v[0])*dr[0] + (C->Flow_data.Grads.grad_v[1])*dr[1] + (C->Flow_data.Grads.grad_v[2])*dr[2]);
      else if(var == 'W')
	*q = (C->Flow_data.State_vec.RhoW)/(C->Flow_data.State_vec.Rhof) + (w_lim)*
	  ((C->Flow_data.Grads.grad_w[0])*dr[0] + (C->Flow_data.Grads.grad_w[1])*dr[1] + (C->Flow_data.Grads.grad_w[2])*dr[2]);
      else if(var == 'E')
	*q = (C->Flow_data.State_vec.RhoE)/(C->Flow_data.State_vec.Rhof) + (E_lim)*
	  ((C->Flow_data.Grads.grad_E[0])*dr[0] + (C->Flow_data.Grads.grad_E[1])*dr[1] + (C->Flow_data.Grads.grad_E[2])*dr[2]);
    }
  }
#else
  /*So we don't want reconstruction*/

  if(var == EOF) {
    sv.Rhof = C -> Flow_data.State_vec.Rhof;
    sv.Pres = C -> Flow_data.State_vec.Pres; 
    sv.RhoU = C -> Flow_data.State_vec.RhoU;
    sv.RhoV = C -> Flow_data.State_vec.RhoV;
    sv.RhoW = C -> Flow_data.State_vec.RhoW;
    sv.RhoE = C -> Flow_data.State_vec.RhoE;
    sv.rhof_prod = C -> Flow_data.State_vec.rhof_prod;
#if AFTERBURN
    sv.rhof_prod_ub = C -> Flow_data.State_vec.rhof_prod_ub;
#endif    
  }
  else if(var == 'P') 
    *q = C->Flow_data.State_vec.Pres; 
  else if(var == 'R')
    *q = C->Flow_data.State_vec.Rhof;
  else if(var == 'U')
    *q = (C->Flow_data.State_vec.RhoU)/(C->Flow_data.State_vec.Rhof);
  else if(var == 'V')
    *q = (C->Flow_data.State_vec.RhoV)/(C->Flow_data.State_vec.Rhof);
  else if(var == 'W')
    *q = (C->Flow_data.State_vec.RhoW)/(C->Flow_data.State_vec.Rhof);
  else if(var == 'E')
    *q = (C->Flow_data.State_vec.RhoE)/(C->Flow_data.State_vec.Rhof);
#endif

  return(sv);
}

/*------------------------------------------------------------------*/
