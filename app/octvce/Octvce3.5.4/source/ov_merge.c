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

/**\file Source file for cell merging functions for cartesian cells*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ov_kernel.h"
#include "ov_merge.h"
#include "ov_lists.h"

#define NOTINCLUDE 0

extern List_leaf Leaves;
extern List_leaf Leaves_tail;

/*------------------------------------------------------------------*/

/**\brief Run through list of all cells and merge them if necessary*/

short int merge_all_cells(void)
{
  List_leaf L = Leaves;

  while(L != NULL)
    {
      if(
#if DETON
	 (L -> cell_loc -> un_det != TRUE) && /*Don't let undetonated cells be merged with other cells*/
#endif
	 (L -> cell_loc -> is_small_cell == TRUE) && (L -> cell_loc -> Merge_data == NULL)) /*Is small cell not already merged*/
	{
	  if(merge_cell(L -> cell_loc) == ERROR)
	    {
	      printf("merge_all_cells(): Error merging cells\n");
	      return(ERROR);
	    }
	}

      L = L -> next;
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Run through list of all cells and unmerge them if necessary*/

void unmerge_all_cells(void)
{
  List_leaf L = Leaves;

  while(L != NULL)
    {
      if((L -> cell_loc -> Merge_data != NULL) && 
	 ((L -> cell_loc -> adapt_flag == REFINE) || (L -> cell_loc -> adapt_flag == REALLY_COARSEN)))
	unmerge_cell(L -> cell_loc, FALSE);
		
      L = L -> next;
    }
}

/*------------------------------------------------------------------*/

/**\brief Recursively find a 'normal' (non-small) cell to merge with - a cell can merge directly
   with a non-merged normal cell, an already merged normal cell or an already merged small cell (this
   is also the order of preference)*/

short int merge_cell(Cart_cell C)
{
  short int i, j, k, merged_dir, merged_quad, opp_dir;
  Cart_cell normal_cell = NULL;
  Cart_cell merged_normal_cell = NULL;
  j = 0; merged_dir = 0; merged_quad = 0; opp_dir = 0;
  /*Now search directions for a normal cell to merge with*/
  for(i = 0; i < NUM_FACES; i++)
    {
      if((C -> face_neighbours[i][0] != NULL) && (C -> face_neighbours[i][1] == NULL) &&
	 ((C -> flux_area[i][0] > 0) || (C -> flux_area[i][1] > 0) || 
	  (C -> flux_area[i][2] > 0) || (C -> flux_area[i][3] > 0))) /*1 neighbour can merge with here*/
	{
	  if(C -> face_neighbours[i][0] -> is_small_cell == FALSE)	      
	    { /*Found a 'normal' cell*/

	      if(C -> face_neighbours[i][0] -> Merge_data == NULL) /*Unmerged normal cell*/
		{
		  normal_cell = C -> face_neighbours[i][0];
		  merged_dir = i; merged_quad = 0;

		  break; 
		}
	      else merged_normal_cell = C -> face_neighbours[i][0]; /*Remember a merged normal cell - may have to merge with it*/
	    }
	}
      else if(C -> face_neighbours[i][1] != NULL) /*4 cells in this direction*/
	{
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
	    {
	      if((C -> flux_area[i][j] > 0) && (C -> face_neighbours[i][j] -> is_small_cell == FALSE))
		{
		  if(C -> face_neighbours[i][j] -> Merge_data == NULL) 
		    {
		      normal_cell = C -> face_neighbours[i][j];
		      merged_dir = i; merged_quad = j;
		      break; 
		    }
		  else merged_normal_cell = C -> face_neighbours[i][j];
		}
	    }

	  if(normal_cell != NULL)
	    break;
	}
    } /*Finished cycling through all directions*/

  if(normal_cell != NULL) /*Most commonly, an unmerged normal cell can be found to merge with*/
    {
      C -> Merge_data = malloc(sizeof(struct merge));

      if(C -> Merge_data == NULL)
	return(ERROR);
      else
	{
	  C -> adjacent_to_normal_cell = TRUE; /*Normal cell is one of this cell's neighbours*/

	  normal_cell -> Merge_data = C -> Merge_data; /*All cells point to same structure*/
	  C -> Merge_data -> Linked_cells = NULL;
	  
	  if(add_to_merge_list_front(normal_cell, &(normal_cell->Merge_data->Linked_cells)) == ERROR) /*Normal cell be list's tail*/
	    return(ERROR);
	  else if(add_to_merge_list_front(C, &(C -> Merge_data -> Linked_cells)) == ERROR)
	    return(ERROR);

	  C -> Merge_data -> cell_volume = (C -> cell_volume) + (normal_cell -> cell_volume);

#if DETON
	  normal_cell -> un_det = FALSE; /*If merge with an undetonated normal cell, let it be detonated*/
#endif

	  /*Finally deallocate any flux vector that exists between these 2 cells*/

	  if(C -> face_neighbours[merged_dir][1] == NULL) /*Normal cell is only cell in direction merged_dir*/
	    {
	      if(C -> Flow_data.Face_fluxes[merged_dir][0] != NULL) /*Flux vector exists between these cells*/
		{
		  switch(merged_dir)
		    {
		    case EST:
		      opp_dir = WST; break;
		    case WST:
		      opp_dir = EST; break;
		    case NTH:
		      opp_dir = STH; break;
		    case STH:
		      opp_dir = NTH; break;
		    case UPR:
		      opp_dir = LWR; break;
		    case LWR:
		      opp_dir = UPR; break;
		    }

		  if((normal_cell -> cell_level) == (C -> cell_level))
		    normal_cell -> Flow_data.Face_fluxes[opp_dir][0] = NULL;
		  else /*Only other option is 1 level coarser*/
		    {
		      /*Need to find corresponding quadrant of neighbour*/

		      switch(merged_dir)
			{
			case NTH:
			  j = ((C -> child_num)-1)/2; break;
			  /*jth quadrant of parent cell that child resides on; normal_cell's quadrant also jth quadrant*/
			case STH:
			  j = (C -> child_num)/2; break;
			case EST:
			  if(C -> child_num <= 3) /*Cells 2 and 3*/
			    j = (C -> child_num)-2;
			  else j = (C -> child_num)-4; /*Cells 6 and 7*/
			  break;
			case WST:
			  if(C -> child_num <= 1) /*Cells 0 and 1*/
			    j = (C -> child_num);
			  else j = (C -> child_num)-2;
			  break;
			case LWR:
			  j = (C -> child_num); break;
			case UPR:
			  j = (C -> child_num)-4; break;
			}
		      
		      normal_cell -> Flow_data.Face_fluxes[opp_dir][j] = NULL; /*normal_cell previously pointed to this
										 flux vector*/
		    }

		  free(C -> Flow_data.Face_fluxes[merged_dir][0]);
		  C -> Flow_data.Face_fluxes[merged_dir][0] = NULL;		  
		} 
	    }
	  else /*Merged with one of 4 cells on this face*/
	    {
	      if(C -> Flow_data.Face_fluxes[merged_dir][merged_quad] != NULL)
		{
		  switch(merged_dir)
		    {
		    case EST:
		      opp_dir = WST; break;
		    case WST:
		      opp_dir = EST; break;
		    case NTH:
		      opp_dir = STH; break;
		    case STH:
		      opp_dir = NTH; break;
		    case UPR:
		      opp_dir = LWR; break;
		    case LWR:
		      opp_dir = UPR; break;
		    }

		  normal_cell -> Flow_data.Face_fluxes[opp_dir][0] = NULL; /*Normal cell is in direction merged_dir, quadrant merged_quad*/
		  free(C -> Flow_data.Face_fluxes[merged_dir][merged_quad]);
		  C -> Flow_data.Face_fluxes[merged_dir][merged_quad] = NULL;
		}
	    }
	}
    }
  else if(merged_normal_cell != NULL) /*At least can merge with already merged normal cell - less common*/
    { 
      C -> adjacent_to_normal_cell = TRUE;

      C -> Merge_data = merged_normal_cell -> Merge_data;

      if(add_to_merge_list_front(C, &(C -> Merge_data -> Linked_cells)) == ERROR)
	return(ERROR);

      C -> Merge_data -> cell_volume = (C -> Merge_data -> cell_volume) + (C -> cell_volume);

#if DETON
      merged_normal_cell -> un_det = FALSE; /*If merge with an undetonated normal cell, let it be detonated*/
#endif

      /*This cell may already belong to a cluster of cells all merged together - deallocate any flux vectors
	shared by this cell with others of this cluster*/

      for(i = 0; i < NUM_FACES; i++)
	{
	  switch(i)
	    {
	    case EST:
	      opp_dir = WST; break;
	    case WST:
	      opp_dir = EST; break;
	    case STH:
	      opp_dir = NTH; break;
	    case NTH:
	      opp_dir = STH; break;
	    case LWR:
	      opp_dir = UPR; break;
	    case UPR:
	      opp_dir = LWR; break;
	    }

	  if((C -> face_neighbours[i][0] != NULL) && (C -> face_neighbours[i][1] == NULL) && 
	     (C -> face_neighbours[i][0] -> Merge_data == C -> Merge_data) && (C -> Flow_data.Face_fluxes[i][0] != NULL))
	    { /*1 neighbour can flux with in this direction but also belongs to cluster of merged cells - deallocate flux vector*/
	      
	      if((C -> face_neighbours[i][0] -> cell_level) == (C -> cell_level))
		C -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][0] = NULL;
	      else
		{
		  switch(i)
		    {
		    case NTH:
		      j = ((C -> child_num)-1)/2; break;
		    case STH:
		      j = (C -> child_num)/2; break;
		    case EST:
		      if(C -> child_num <= 3)
			j = (C -> child_num)-2;
		      else j = (C -> child_num)-4;
		      break;
		    case WST:
		      if(C -> child_num <= 1) 
			j = (C -> child_num);
		      else j = (C -> child_num)-2;
		      break;
		    case LWR:
		      j = (C -> child_num); break;
		    case UPR:
		      j = (C -> child_num)-4; break;
		    }
		  
		  C -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][j] = NULL; 
		}

	      free(C -> Flow_data.Face_fluxes[i][0]);
	      C -> Flow_data.Face_fluxes[i][0] = NULL;
	    }	  
	  else if(C -> face_neighbours[i][1] != NULL)
	    {
	      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		{
		  if((C -> face_neighbours[i][j] -> Merge_data == C -> Merge_data) && (C -> Flow_data.Face_fluxes[i][j] != NULL))
		    {				  
		      C -> face_neighbours[i][j] -> Flow_data.Face_fluxes[opp_dir][0] = NULL;
		      free(C -> Flow_data.Face_fluxes[i][j]);
		      C -> Flow_data.Face_fluxes[i][j] = NULL;
		    }
		}
	    }
	}
    }
  else /*Very rare situation when all neighbours are small/solid/non-existent e.g. corner small cell - merge with
	 a cluster of linked cells*/
    { 
      Cart_cell small_cell = NULL; /*Find another small cell to join with that should be (possibly) joined to a group of small
				     cells merged with a normal cell*/
      
      for(i = 0; i < NUM_FACES; i++)
	{
	  if((C -> face_neighbours[i][0] != NULL) && (C -> face_neighbours[i][1] == NULL) && 
	     (C -> face_neighbours[i][0] -> cell_marked != TRUE) &&
	     ((C -> flux_area[i][0] > 0) || (C -> flux_area[i][1] > 0) || 
	      (C -> flux_area[i][2] > 0) || (C -> flux_area[i][3] > 0)))
	    {
	      small_cell = C -> face_neighbours[i][0];
	      break;
	    }
	  else if(C -> face_neighbours[i][1] != NULL)
	    {
	      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		{
		  if((C -> face_neighbours[i][j] -> cell_marked != TRUE) && (C -> flux_area[i][j] > 0))
		    {
		      small_cell = C -> face_neighbours[i][j];
		      break;
		    }
		}

	      if(small_cell != NULL)
		break;
	    }
	}

      C -> cell_marked = TRUE; /*To prevent another small cell neighbour to this cell returning to this cell*/

      if((small_cell != NULL) && (small_cell -> Merge_data == NULL))
	{
	  if(merge_cell(small_cell) == ERROR) /*Recursively find a normal cell to merge with*/
	    return(ERROR);
	}
      
      if((small_cell != NULL) && (small_cell -> Merge_data != NULL)) /*A normal cell to merge with successfully found*/
	{
	  C -> cell_marked = UNKNOWN; /*Reset*/
	  
	  C -> adjacent_to_normal_cell = UNKNOWN; /*Reset - normal cell not any of this cell's neighbours*/

	  C -> Merge_data = small_cell -> Merge_data;

	  if(add_to_merge_list_front(C, &(C -> Merge_data -> Linked_cells)) == ERROR)
	    return(ERROR);

	  C -> Merge_data -> cell_volume = (C -> Merge_data -> cell_volume) + (C -> cell_volume);

#if DETON
	  small_cell -> un_det = TRUE;
#endif

	  /*Now finally deallocate flux vectors if required*/

	  for(i = 0; i < NUM_FACES; i++)
	    {
	      switch(i)
		{
		case EST:
		  opp_dir = WST; break;
		case WST:
		  opp_dir = EST; break;
		case STH:
		  opp_dir = NTH; break;
		case NTH:
		  opp_dir = STH; break;
		case LWR:
		  opp_dir = UPR; break;
		case UPR:
		  opp_dir = LWR; break;
		}
	      
	      if((C -> face_neighbours[i][0] != NULL) && (C -> face_neighbours[i][1] == NULL) && 
		 (C -> face_neighbours[i][0] -> Merge_data == C -> Merge_data) && (C -> Flow_data.Face_fluxes[i][0] != NULL))
		{		  
		  if((C -> face_neighbours[i][0] -> cell_level) == (C -> cell_level))
		    C -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][0] = NULL;
		  else
		    {
		      switch(i)
			{
			case NTH:
			  j = ((C -> child_num)-1)/2; break;
			case STH:
			  j = (C -> child_num)/2; break;
			case EST:
			  if(C -> child_num <= 3)
			    j = (C -> child_num)-2;
			  else j = (C -> child_num)-4;
			  break;
			case WST:
			  if(C -> child_num <= 1) 
			    j = (C -> child_num);
			  else j = (C -> child_num)-2;
			  break;
			case LWR:
			  j = (C -> child_num); break;
			case UPR:
			  j = (C -> child_num)-4; break;
			}
		      
		      C -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][j] = NULL; 
		    }
		  
		  free(C -> Flow_data.Face_fluxes[i][0]);
		  C -> Flow_data.Face_fluxes[i][0] = NULL;		  
		}	  
	      else if(C -> face_neighbours[i][1] != NULL)
		{
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		    {
		      if((C -> face_neighbours[i][j] -> Merge_data == C -> Merge_data) && (C -> Flow_data.Face_fluxes[i][j] != NULL))
			{				  
			  C -> face_neighbours[i][j] -> Flow_data.Face_fluxes[opp_dir][0] = NULL;
			  free(C -> Flow_data.Face_fluxes[i][j]);
			  C -> Flow_data.Face_fluxes[i][j] = NULL;
			}
		    }
		}
	    }	  
	}
      else 
	{ /*This should only happen in the rarest of situations when a small cell or group of small cells is isolated
	    from the rest of the flow - make the cell a solid one*/
	  
	  for(i = 0; i < NUM_FACES; i++)
	    {
	      switch(i)
		{
		case EST:
		  opp_dir = WST; break;
		case WST:
		  opp_dir = EST; break;
		case STH:
		  opp_dir = NTH; break;
		case NTH:
		  opp_dir = STH; break;
		case LWR:
		  opp_dir = UPR; break;
		case UPR:
		  opp_dir = LWR; break;
		}
	      
	      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		C -> flux_area[i][j] = 0;

	      if(C -> face_neighbours[i][0] != NULL)
		{
		  if(C -> face_neighbours[i][0] -> cell_level == (C -> cell_level) - 1) /*Coarser neighbour*/
		    {
		      switch(i)
			{
			case NTH:
			  j = ((C -> child_num)-1)/2; break;
			case STH:
			  j = (C -> child_num)/2; break;
			case EST:
			  if(C -> child_num <= 3)
			    j = (C -> child_num)-2;
			  else j = (C -> child_num)-4;
			  break;
			case WST:
			  if(C -> child_num <= 1) 
			    j = (C -> child_num);
			  else j = (C -> child_num)-2;
			  break;
			case LWR:
			  j = (C -> child_num); break;
			case UPR:
			  j = (C -> child_num)-4; break;
			}
		      
		      C -> face_neighbours[i][0] -> flux_area[opp_dir][j] = 0;
		    }
		  else if(C -> face_neighbours[i][0] -> cell_level == C -> cell_level) /*Same level neighbour*/
		    {
		      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
			C -> face_neighbours[i][0] -> flux_area[opp_dir][j] = 0;
		    }
		  else if(C -> face_neighbours[i][0] -> cell_level == (C -> cell_level) + 1) /*Finer neighbour*/
		    {
		      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Cycle through each cell on each face*/
			{
			  for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
			    C -> face_neighbours[i][j] -> flux_area[opp_dir][k] = 0;
			}
		    }
		}
	    }

	  C -> cell_type = SOLID;
	  C -> cell_volume = 0;
	  C -> cell_marked = UNKNOWN;
	  delete_from_leaf_list(C, &Leaves, &Leaves_tail); 
	} 
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Unmerge a cell from a cluster or linked cells, and possibly destroy this cluster (unlink all cells) if
   this cell is the 'normal' cell of the cluster or there's only 1 cell left in it.  Will also recursively check if
   other (non-adjacent to normal cell) small cells part of the same cluster need to be unmerged too - will unmerge all
   cells that need unmerging in most realistic circumstances*/

short int unmerge_cell(Cart_cell C, short int check_if_connected)
{  
  if(C -> is_small_cell == FALSE) /*Is the 'normal' cell of this merged cluster - destroy the merged cell (unlink all cells)*/
    {
      List_merge mlist, temp;
      Merge mdata;
      
      mlist = C -> Merge_data -> Linked_cells; /*Head of list*/
      temp = mlist;
      mdata = C -> Merge_data;
      
      while(temp != NULL) /*Destroy the list first*/
	{
	  temp = temp -> next;
	  
	  mlist -> cell_loc -> Merge_data = NULL;
	  mlist -> cell_loc -> Merge_list_loc = NULL;
	  mlist -> cell_loc -> adjacent_to_normal_cell = UNKNOWN;
	  
	  free(mlist);
	  
	  mlist = temp;
	}
      
      free(mdata); /*Now destroy the structure containing the list*/

      return(NOERROR);
    }
  else /*Is a small cell on the merged cluster - it may not be necessary to destroy whole merged cluster*/
    {
      short int i, j;
      short int is_connected = FALSE;
      C -> cell_marked = TRUE; /*Currently 'investigating' this cell*/
      
      /*Check if this cell adjacent to any other small cells in same cluster - if as a result of unmerging this
	cell other small cells no longer are connected to the cluster, we must also unmerge them first.
	
	In certain extreme situations some cells will be unmerged that won't have to, but I feel this is the
	most cost-effective approach that accomplishes the most legitimate unmerging.  Unnecessarily unmerged cells
	can always be re-merged or merged with other normal cells*/

      for(i = 0; i < NUM_FACES; i++)
	{
	  if((C -> face_neighbours[i][0] != NULL) && (C -> face_neighbours[i][1] == NULL) && 
	     (C -> face_neighbours[i][0] -> Merge_data == C -> Merge_data) && (C -> face_neighbours[i][0] -> is_small_cell == TRUE) &&
	     (C -> face_neighbours[i][0] -> cell_marked != TRUE) &&
	     ((C -> flux_area[i][0] > 0) || (C -> flux_area[i][1] > 0) || 
	      (C -> flux_area[i][2] > 0) || (C -> flux_area[i][3] > 0)))
	    { /*See if any small neighbour is adjacent to the normal cell, if so this cell can definitely
		stay merged if check_if_connected TRUE*/
	      
	      if(C -> face_neighbours[i][0] -> adjacent_to_normal_cell == TRUE)
		is_connected = TRUE;
	      else
		{
		  if(unmerge_cell(C -> face_neighbours[i][0], TRUE) == CONNECTED)
		    is_connected = TRUE; /*This cell is adjacent to another small cell connected to the normal cell*/
		}
	      
	      if((check_if_connected == TRUE) && (is_connected == TRUE))
		break; /*Found that the cell is connected (it can be spared unmerging if check_if_connected FALSE; if 
			 check_if_connected FALSE we want to check all directions)*/
	    }
	  else if(C -> face_neighbours[i][1] != NULL)
	    {
	      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		{
		  if((C->face_neighbours[i][j]->Merge_data == C->Merge_data) && (C->face_neighbours[i][j]->is_small_cell == TRUE) &&
		     (C -> face_neighbours[i][j] -> cell_marked != TRUE) && (C -> flux_area[i][j] > 0))
		    {
		      if(C -> face_neighbours[i][j] -> adjacent_to_normal_cell == TRUE)
			is_connected = TRUE;
		      else
			{
			  if(unmerge_cell(C -> face_neighbours[i][j], TRUE) == CONNECTED)
			    is_connected = TRUE;
			}
		      
		      if((check_if_connected == TRUE) && (is_connected == TRUE))
			break;
		    }
		}
	      
	      if((check_if_connected == TRUE) && (is_connected == TRUE))
		break;
	    }
	}
      
      C -> cell_marked = UNKNOWN; /*Reset*/
      
      if((is_connected == FALSE) || (check_if_connected == FALSE)) /*Must unmerge this cell*/
	{
	  delete_from_merge_list(C, &(C -> Merge_data -> Linked_cells));
	  
	  if((C -> Merge_data -> Linked_cells -> prev == NULL) && (C -> Merge_data -> Linked_cells -> next == NULL))
	    { /*1 node left on list (will be the normal cell) - destroy the cluster*/

	      C -> Merge_data -> Linked_cells -> cell_loc -> Merge_list_loc = NULL;
	      C -> Merge_data -> Linked_cells -> cell_loc -> Merge_data = NULL;
	      C -> Merge_data -> Linked_cells -> cell_loc -> adjacent_to_normal_cell = UNKNOWN; /*Actually a redundant step*/
	      
	      free(C -> Merge_data -> Linked_cells);
	      free(C -> Merge_data);
	    }
	  
	  C -> Merge_data = NULL;
	  C -> adjacent_to_normal_cell = UNKNOWN;	  
	}
      else if(is_connected == TRUE) /*No need to unmerge this cell - it's connected to the normal cell*/
	return(CONNECTED);
    }

  return(NOT_CONNECTED);  
}

/*------------------------------------------------------------------*/
