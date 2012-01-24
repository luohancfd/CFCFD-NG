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

/**\file Source file for establishing connectivities between cells - flux, neighbour and vertex connectivities*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ov_kernel.h"
#include "ov_connect.h"
#include "ov_lists.h"
#include "ov_setgeom.h"

#define NOTINCLUDE 0

#if !NOTINCLUDE
int atime = 0;
extern int stepnow;
#endif

extern short int min_refinement_level;
extern short int refining_root_now;
extern double root_scale;
extern List_vtx_glob Vtxlist;
extern List_vtx_glob Vtxlist_tail;
extern List_vtx_glob *Vtx_heads;
extern List_vtx_glob *Vtx_tails;
extern List_leaf Leaves;
extern short int adapt_in_parallel;
extern int number_of_threads;

/*------------------------------------------------------------------*/

/**\brief Establish neighbouring connectivities in serial after a cell is refined once*/

void establish_neighbour_after_refine(Cart_cell child_cells[], short int directions) 
{                                            						
  short int i, j, k, y, z, dir, opp_dir; 
  Cart_cell Parent = child_cells[0] -> parent; /*Parent for all children is the same*/
 
#if NOTINCLUDE
 if((Parent->cell_level==5) && (Parent->centroid[0]==0.015625) && (Parent->centroid[1]==-0.015625) && 
	(Parent->centroid[2] == -0.203125)) {
 	printf("We do want to currently do this cell thanks\n");
	if(directions == ALL)
	 printf("directions is ALL\n");
 }
#endif
 
  short int dirs[3];
  short int num_faces;
  y = 0; z = 0; opp_dir = 0; num_faces = 0;
  if(directions == ALL)
    num_faces = NUM_FACES;
  else if(directions == POSITIVE)
    {
      dirs[0] = EST; dirs[1] = NTH; dirs[2] = UPR;
      num_faces = 3;
    }
  else if(directions == NEGATIVE)
    {
      dirs[0] = WST; dirs[1] = STH; dirs[2] = LWR;
      num_faces = 3;
    }

  if(directions != NEGATIVE) /*Can establish internal connectivites fairly easily*/
    {
      for(i = 0; i < NUM_FACES; i++)
	{	  
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
	    {
	      switch(i)
		{
		case NTH:
		  y = 2*j+1; z = y-1; break; /*z is the 'opposite' direction from y*/
		case STH:
		  y = 2*j; z = y+1; break;
		case EST:
		  if(j <= 1)
		    y = j+2;
		  else y = j+4;
		  z = y-2; break;
		case WST:
		  if(j <= 1)
		    y = j;
		  else y = j+2;
		  z = y+2; break;
		case LWR: 
		  y = j; z = y+4; break;
		case UPR: 
		  y = j+4; z = y-4; break;
		}

	      child_cells[z] -> face_neighbours[i][0] = child_cells[y];
	      /*Child cells on ith face (cells y) are in direction i with respect to child cells on the opposite face (cells z)*/
	    }
	}
    }

  /*Now go through parent's exterior faces and update connectivities*/

  for(dir = 0; dir < num_faces; dir++) /*Cycle through each face/direction*/
    {
      if(num_faces == 3)
	i = dirs[dir];
      else i = dir;

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
      
      /*Now go and update children connectivities*/

      if(Parent -> face_neighbours[i][0] == NULL) /*Parent has no neighbour in direction i*/
	{
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Cycle through child cells on ith face*/
	    {
	      switch(i)
		{
		case NTH: /*Cells on north face are 1,3,5,7*/
		  y = 2*j+1; break;
		case STH: 
		  y = 2*j; break;
		case EST:
		  if(j <= 1)
		    y = j+2;
		  else y = j+4;
		  break;
		case WST:
		  if(j <= 1)
		    y = j;
		  else y = j+2;
		  break;
		case LWR:
		  y = j; break;
		case UPR:
		  y = j+4; break;
		}
	      
	      child_cells[y] -> face_neighbours[i][0] = NULL;
	    }
	}
      else if (Parent -> face_neighbours[i][1] == NULL)
	{ 
	  /*Only 1 neighbour in this direction*/

	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Cycle through child cells on ith face*/
	    {
	      switch(i)
		{
		case NTH: 
		  y = 2*j+1; break;
		case STH: 
		  y = 2*j; break;
		case EST:
		  if(j <= 1)
		    y = j+2;
		  else y = j+4;
		  break;
		case WST:
		  if(j <= 1)
		    y = j;
		  else y = j+2;
		  break;
		case LWR:
		  y = j; break;
		case UPR:
		  y = j+4; break;
		}

	      child_cells[y] -> face_neighbours[i][0] = Parent -> face_neighbours[i][0];

	      if(Parent -> face_neighbours[i][0] -> cell_level == Parent -> cell_level)
		Parent -> face_neighbours[i][0] -> face_neighbours[opp_dir][j] = child_cells[y];
	      
	      /*Parent's only neighbour in this direction updates its connectivities with the parent's children*/
	    }
	} 
      else /*So parent has more than 1 neighbour in this direction*/
	{
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Cycle through parent's children on this face*/
	    {
	      switch(i)
		{
		case NTH: 
		  y = 2*j+1; break;
		case STH: 
		  y = 2*j; break;
		case EST:
		  if(j <= 1)
		    y = j+2;
		  else y = j+4;
		  break;
		case WST:
		  if(j <= 1)
		    y = j;
		  else y = j+2;
		  break;
		case LWR:
		  y = j; break;
		case UPR:
		  y = j+4; break;
		}
	      
	      if(Parent -> face_neighbours[i][j] -> children[0] == NULL) /*Neighbour doesn't have children*/
		{
		  child_cells[y] -> face_neighbours[i][0] = Parent -> face_neighbours[i][j];
		  /*Child cell updates with its corresponding neighbour*/
		  
		  Parent -> face_neighbours[i][j] -> face_neighbours[opp_dir][0] = child_cells[y];
		  /*Parent's 4 neighbours update their connectivities with children*/
		}
	      else 
		{ /*Neighbour actually has children of its own; the child cell should update with these children*/

		  for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
		    {
		      switch(i)
			{
			case NTH:
			  z = 2*k; break; /*This is the opposite direction to north*/
			case STH:
			  z = 2*k+1; break;
			case EST:
			  if(k <= 1)
			    z = k;
			  else z = k+2;
			  break;
			case WST:
			  if(k <= 1)
			    z = k+2;
			  else z = k+4;
			  break;
			case LWR:
			  z = k+4; break;
			case UPR:
			  z = k; break;
			}

		      child_cells[y] -> face_neighbours[i][k] = Parent -> face_neighbours[i][j] -> children[z];
		      /*The child cell updates with the opposite-faced children of the neighbour*/

		      child_cells[y] -> face_neighbours[i][k] -> face_neighbours[opp_dir][0] = child_cells[y];
		      /*The neighbour's children update their connectivities with the child*/

		      if(child_cells[y] -> face_neighbours[i][k] -> children[0] != NULL)
			assign_new_face_neighbour(child_cells[y] -> face_neighbours[i][k], opp_dir, child_cells[y]);
		      /*Recursively traverse all cells adjacent to this face*/
		    }

		  /*Neighbour's parent updates with child cell too*/		  
		  child_cells[y] -> face_neighbours[i][0] -> parent -> face_neighbours[opp_dir][0] = child_cells[y];
		}	      
	    }
	}
    }

  /*All neighbouring connectivities finally established*/
}

/*------------------------------------------------------------------*/

/**\brief Establish neighbours after coarsening a cell*/

void establish_neighbour_after_coarsen(Cart_cell C, short int directions)
{
  short int i, j, k, dir, opp_dir;
  
  short int dirs[3];
  short int num_faces;
  k = 0; opp_dir = 0; num_faces = 0;
  if(directions == ALL)
    num_faces = NUM_FACES;
  else if(directions == POSITIVE)
    {
      dirs[0] = EST; dirs[1] = NTH; dirs[2] = UPR;
      num_faces = 3;
    }
  else if(directions == NEGATIVE)
    {
      dirs[0] = WST; dirs[1] = STH; dirs[2] = LWR;
      num_faces = 3;
    }

  for(dir = 0; dir < num_faces; dir++) 
    {
      if(num_faces == 3)
	i = dirs[dir];
      else i = dir;

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
			      
      /*First test if children in this direction have a coarser neighbour*/
			    
      if((i == EST) || (i == NTH) || (i == UPR))
	j = 7; /*Can look at child 7*/
      else j = 0;
  			
      if(C -> children[j] -> face_neighbours[i][0] == NULL)
	{
	  for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
	    C -> face_neighbours[i][k] = NULL;
	}
      else if(C -> children[j] -> face_neighbours[i][0] -> cell_level <= C -> cell_level)
	{ /*Parent has only 1 neighbour (or none) in this direction*/
			  
	  C -> face_neighbours[i][0] = C -> children[j] -> face_neighbours[i][0];
	  /*Parent's face neighbour is its children's face neighbour*/
				  
	  for(k = 1; k < MAX_NUM_NEIGHBOURS; k++) /*Parent only stores 1 face neighbour*/ 	  
	    C -> face_neighbours[i][k] = NULL;
				  
	  /*Must update the corresponding neighbour*/
	  	  
	  if(C -> face_neighbours[i][0] -> cell_level == C -> cell_level)
	    { /*Neighbour is at same level; can update connectivity, else don't*/
					  
	      C -> face_neighbours[i][0] -> face_neighbours[opp_dir][0] = C;
					      
	      for(k = 1; k < MAX_NUM_NEIGHBOURS; k++)
		C -> face_neighbours[i][0] -> face_neighbours[opp_dir][k] = NULL;
	    }
	}
      else /*So each child might have 1 neighbour each, or more than 1 neighbour per child - but only
	     4 neighbours are stored still*/
	{	
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
	    {
	      switch(i)
		{
		case NTH: 
		  k = 2*j+1; break;
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
		  else k = j+2;
		  break;
		case LWR: 
		  k = j; break;
		case UPR: 
		  k = j+4; break;
		}
				      
	      if(C -> children[k] -> face_neighbours[i][1] == NULL) 
		{ /*This child only has 1 neighbour*/

		  /*Set cell's face neighbour in this quadrant*/
		  C -> face_neighbours[i][j] = C -> children[k] -> face_neighbours[i][0];
					  
		  /*Update corresponding neighbour*/
		  C -> face_neighbours[i][j] -> face_neighbours[opp_dir][0] = C;
		}
	      else /*This child stores more than 1 neighbour already*/
		{
		  /*Set the cell's face neighbour in this quadrant to be the child's neighbour's parent*/
		  C -> face_neighbours[i][j] = C -> children[k] -> face_neighbours[i][0] -> parent;
					  
		  /*Update the corresponding neighbour's parent*/
		  C -> face_neighbours[i][j] -> face_neighbours[opp_dir][0] = C;
					  
		  /*Now finally let all descendants of the neighbour on this face point to parent*/
		  assign_new_face_neighbour(C -> face_neighbours[i][j], opp_dir, C);
		}
	    }					
	}				  
    } /*Finished updating neighbouring connectivities*/
}


/*------------------------------------------------------------------*/

/**\brief Check a cell's non-sibling neighbours (before it's refined) to see if they need to be refined first to ensure grid smoothness*/

short int check_neighbour_for_refine(Cart_cell C, Cart_cell neighbours[]) 
{ 
  short int i, dir;
  short int n = 0; /*No. neighbours that need refinement*/

  if(C -> parent != NULL) 
    {
      /*Determine what child of parent it is - will affect what neighbours can be refined*/

      short int child_num = C -> child_num;

      /*Neighbours that need refinement are those neighbouring cells which are already a coarser level than this cell's and can
	flux with this cell (only non-sibling neighbours as sibling neighbours won't be coarser than cell).*/	
      
      short int non_sibling_neighbours[3]; /*At most 3 non-sibling neighbours*/

      switch(child_num)
	{
	case 0: 
	  non_sibling_neighbours[0] = STH; non_sibling_neighbours[1] = WST; non_sibling_neighbours[2] = LWR;
	  /*Directions of non-sibling neighbours of cell 0*/
	  break;
	case 1: 
	  non_sibling_neighbours[0] = NTH; non_sibling_neighbours[1] = WST; non_sibling_neighbours[2] = LWR; break;
	case 2: 
	  non_sibling_neighbours[0] = STH; non_sibling_neighbours[1] = EST; non_sibling_neighbours[2] = LWR; break;
	case 3: 
	  non_sibling_neighbours[0] = NTH; non_sibling_neighbours[1] = EST; non_sibling_neighbours[2] = LWR; break;
	case 4: 
	  non_sibling_neighbours[0] = STH; non_sibling_neighbours[1] = WST; non_sibling_neighbours[2] = UPR; break;
	case 5: 
	  non_sibling_neighbours[0] = NTH; non_sibling_neighbours[1] = WST; non_sibling_neighbours[2] = UPR; break;
	case 6: 
	  non_sibling_neighbours[0] = STH; non_sibling_neighbours[1] = EST; non_sibling_neighbours[2] = UPR; break;
	case 7: 
	  non_sibling_neighbours[0] = NTH; non_sibling_neighbours[1] = EST; non_sibling_neighbours[2] = UPR; break;
	}
          
      if((adapt_in_parallel == FALSE) || (refining_root_now == TRUE))
	{  
	  for(i = 2; i >= 0; i--) 
	    {
	      dir = non_sibling_neighbours[i];

	      if((C -> flux_area[dir][0] != 0) || (C -> flux_area[dir][1] != 0) || (C -> flux_area[dir][2] != 0) || 
		 (C -> flux_area[dir][3] != 0))
		{ /*Even if area is unknown (or slightly non-zero), we should still refine anyway*/

		  if((C -> face_neighbours[dir][0] != NULL) && (C -> face_neighbours[dir][0] -> cell_level < C -> cell_level))
		    { /*Non-sibling neighbour's level coarser than cell's*/
		  
		      neighbours[n] = C -> face_neighbours[dir][0];
		      n++;
		    }
		}	  
	    }
	}
      else /*For parallel adaptation we merely want to flag this cell for refinement*/
	{
	  for(i = 2; i >= 0; i--) 
	    {
	      dir = non_sibling_neighbours[i];

	      if((C -> flux_area[dir][0] != 0) || (C -> flux_area[dir][1] != 0) || (C -> flux_area[dir][2] != 0) || 
		 (C -> flux_area[dir][3] != 0))
		{ 
		  if((C -> face_neighbours[dir][0] != NULL) && (C -> face_neighbours[dir][0] -> cell_level < C -> cell_level))
		    { 
		      C -> face_neighbours[dir][0] -> adapt_flag = REFINE;
		      check_neighbour_for_refine(C -> face_neighbours[dir][0], neighbours);
		    }
		}	  
	    }
	}
    }
  return(n);
}

/*------------------------------------------------------------------*/

/**\brief Check to see if a cell can genuinely be coarsened (i.e. children don't have children, geometric data OK, and
   neighbour connectivities allow it).  Returns FALSE if can't be coarsened, TRUE if can.*/

short int check_cell_for_coarsening(Cart_cell C)
{
  short int i, j, k, x, y; 
  k = 0; 
  for(i = 0; i < MAX_NUM_CHILDREN; i++) /*Check if children have children*/
    {
      if(C -> children[i] -> children[0] != NULL)
	{
	  for(y = 0; y < MAX_NUM_CHILDREN; y++) 
	    {
	      if(C -> children[y] -> adapt_flag != REFINE)
		C -> children[y] -> adapt_flag = UNKNOWN; /*Reset*/
	    }

	  return(FALSE); /*Children have children - cannot coarsen cell*/
	}
    }	      
  
  if(set_parent_geom_data(C -> children) == ERROR)
    return(ERROR); /*Must ensure all geometric data known before coarsening*/
  	  
  if((C -> additional_refine_times > 0) || (C -> cell_level < min_refinement_level))
    {
      for(y = 0; y < MAX_NUM_CHILDREN; y++)
	{
	  if(C -> children[y] -> adapt_flag != REFINE)
	    C -> children[y] -> adapt_flag = UNKNOWN;
	}

      return(FALSE); /*Too coarse for 'geometrical' reasons*/
    }
  else /*Can coarsen for geometry reasons*/  
    {	         	    		  
      for(i = 0; i < NUM_FACES; i++) /*Cycle through each face/direction - check neighbour level differences*/
	{
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Cycle through each child on each face*/
	    {
	      if(C -> flux_area[i][j] > 0) /*Only care for fluxable neighbours*/
		{ 
		  switch(i)
		    {
		    case NTH: /*Children 1,3,5,7 on north face*/
		      k = 2*j+1; break;
		    case STH: /*Children 0,2,4,6*/
		      k = 2*j; break;
		    case EST: /*Children 2,3,6,7*/
		      if(j <= 1)
			k = j+2;
		      else k = j+4;
		      break;
		    case WST: /*Children 0,1,4,5*/
		      if(j <= 1)
			k = j;
		      else k = j+2;
		      break;
		    case LWR: /*Children 0,1,2,3*/
		      k = j; break;
		    case UPR: /*Children 4,5,6,7*/
		      k = j+4; break;
		    }

		  if(C -> children[k] -> face_neighbours[i][1] != NULL)
		    { 
		      /*Test if any of the child's neighbours through this quadrant are marked for coarsening*/

		      for(x = 0; x < MAX_NUM_CHILDREN; x++) /*See if the whole neighbour can be coarsened - a necessary condition
							      would be that all non-solid children have COARSEN/REALLY_COARSEN flags*/
			{
			  if((C -> children[k] -> face_neighbours[i][1] -> parent -> children[x] -> cell_type != SOLID) && 
			     (C -> children[k] -> face_neighbours[i][1] -> parent -> children[x] -> adapt_flag != COARSEN) &&
			     (C -> children[k] -> face_neighbours[i][1] -> parent -> children[x] -> adapt_flag != REALLY_COARSEN))
			    {
			      for(y = 0; y < MAX_NUM_CHILDREN; y++)
				if(C -> children[k] -> face_neighbours[i][1] -> parent -> children[y] -> adapt_flag != REFINE)
				  C -> children[k] -> face_neighbours[i][1] -> parent -> children[y] -> adapt_flag = UNKNOWN;
			
			      for(y = 0; y < MAX_NUM_CHILDREN; y++)
				if(C -> children[y] -> adapt_flag != REFINE)
				  C -> children[y] -> adapt_flag = UNKNOWN;
	    			     
			      return(FALSE); /*Couldn't coarsen neighbour - unset flags so we don't check again*/
			    }
			}
			  
		      if(x == MAX_NUM_CHILDREN) /*Now check if the neighbour's parent can be coarsened*/
			{
			  check_cell_for_coarsening(C -> children[k] -> face_neighbours[i][1] -> parent);

			  /*By now we can discern if neighbours of children can be coarsened (they all should have REALLY_COARSEN flags)*/
		      
			  if(C -> children[k] -> face_neighbours[i][1] -> adapt_flag != REALLY_COARSEN)
			    {
			      for(y = 0; y < MAX_NUM_CHILDREN; y++)
				if(C -> children[y] -> adapt_flag != REFINE) 
				  C -> children[y] -> adapt_flag = UNKNOWN;
			  
			      return(FALSE); /*Child sees more than 1 face neighbour in this direction - can't coarsen*/
			    }			  
			}		      								      
		    }

		  /*Perhaps a child's neighbour is the same level as the child and is flagged for refinement?*/

		  if((C -> children[k] -> face_neighbours[i][0] != NULL) && 
		     (C -> children[k] -> face_neighbours[i][0] -> cell_level == C -> children[k] -> cell_level) && 
		     (C -> children[k] -> face_neighbours[i][0] -> adapt_flag == REFINE))
		    {
		      for(y = 0; y < MAX_NUM_CHILDREN; y++)
			{
			  if(C -> children[y] -> adapt_flag != REFINE)
			    C -> children[y] -> adapt_flag = UNKNOWN;
			}

		      return(FALSE);
		    }
		}
	    }
	}
    }

  /*Extra step for parallel coarsening - cell can really be coarsened*/
  for(i = 0; i < MAX_NUM_CHILDREN; i++)
    C -> children[i] -> adapt_flag = REALLY_COARSEN;

  return(TRUE);  
}

/*------------------------------------------------------------------*/

/**\brief Recursively traverses a cell and its descendants and lets the children cells on a specified
 face point to a new cell 'New_neighb' in that direction.  Assumes descendants cells previously only pointed to 1 
 neighbour in this direction to begin with*/

void assign_new_face_neighbour(Cart_cell C, short int dir, Cart_cell New_neighb)
{
  short int i, child;
  child = 0;
  for(i = 0; i < MAX_NUM_NEIGHBOURS; i++) /*Cycle through all children on this face*/
    {
      switch(dir) /*Get the child on the corresponding face*/
	{
	case EST:
	  if(i <= 1)
	    child = i+2;
	  else child = i+4;
	  break;
	case WST:
	  if(i <= 1)
	    child = i;
	  else child = i+2;
	  break;
	case NTH:
	  child = 2*i+1; break;
	case STH:
	  child = 2*i; break;
	case UPR:
	  child = i+4; break;
	case LWR:
	  child = i; break;
	}

      C -> children[child] -> face_neighbours[dir][0] = New_neighb; /*Neighbour has been reassigned*/
      
      if(C -> children[child] -> children[0] != NULL)
	assign_new_face_neighbour(C -> children[child], dir, New_neighb); /*Now do the same thing with the child's children*/
    }
}

/*------------------------------------------------------------------*/

/**\brief Assign verticies in serial to children after refining a cell*/

short int assign_vtx_after_refine(Cart_cell child_cells[], short int phase, List_vtx_glob *Vtxhead, List_vtx_glob *Vtxtail) 
{
  short int i, j, dir, cell_level, child, mid_vertex;
  double cell_length, diag_cen[3], diag_shift[3], diag_bbox[2][3];
  Cart_cell Diag_cell;
  Vertex Vtx;
    
  Cart_cell Parent = child_cells[0] -> parent;
  cell_level = child_cells[0] -> cell_level;
  cell_length = CALC_CELL_EDGE_LENGTH(child_cells[0]);
  child = 0; mid_vertex = 0; Vtx = NULL;
  /*'INITIAL' PHASE*/
  
  if((phase == ALL) || (phase == ALLOC_CHILD))
    {
      /*Assign volume-centre vertex (in the centre of the parent, after child creation)*/
      
      Vtx = malloc(sizeof(struct vertex));
      
      if(Vtx == NULL)
	return(ERROR);
      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

      if(Vtx -> vtx_nums == NULL)
	return(ERROR);

      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR) /*Add to global list of verticies*/
	return(ERROR);
       
      get_vtx_loc(child_cells[0], 7, Vtx -> loc); /*Is vertex 7 of child 0*/

      for(i = 0; i < MAX_NUM_CHILDREN; i++)
	{
	  child_cells[i] -> verticies[7 - i] = Vtx;
	  Vtx -> Leaf_cells[i] = child_cells[i];
	}

      /*First assign parental verticies to respective children - this can be done safely for all corners as cells in each
	vertex's Leaf_cells array occupy a unique position*/

      for(i = 0; i < MAX_NUM_CHILDREN; i++)
	{
	  child_cells[i] -> verticies[i] = Parent -> verticies[i]; /*Corner verticies reassigned*/
	  child_cells[i] -> verticies[i] -> Leaf_cells[7-i] = child_cells[i];

	  /*Note if cell points to vertex i, it will be stored on the 7-i th position on the vertex array*/
	} 
    }
  
  /*'POSITIVE' PHASE*/
  
  if((phase == ALL) || (phase == POSITIVE))
    {
      /*Now assign mid-plane (of parent) verticies to children on 'positive' faces - but see if neighbours already have them*/

      for(i = 0; i < 3; i++) /*Cycle through each face*/
	{
	  dir = 2*i; /*All 'positive' directions*/

	  if(Parent -> face_neighbours[dir][1] != NULL) /*4 neighbours on this face - vertex already exists*/
	    {
	      switch(dir)
		{
		case EST:
		  Vtx = child_cells[7] -> face_neighbours[dir][0] -> verticies[0]; break; /*Works even if more than 1 face neighbour here*/
		case NTH:
		  Vtx = child_cells[7] -> face_neighbours[dir][0] -> verticies[0]; break;
		case UPR:
		  Vtx = child_cells[7] -> face_neighbours[dir][0] -> verticies[0]; break;
		}
	    }
	  else /*Must allocate vertex - it doesn't exist yet*/
	    {
	      Vtx = NULL;
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(j = 0; j < MAX_NUM_CHILDREN; j++)
		Vtx -> Leaf_cells[j] = NULL; /*Initialize to NULL*/

	      /*We can get the vertex locations easily - simply shift centroid of parent to appropriate face*/

	      if(dir == EST) 
		{
		  Vtx -> loc[0] = (Parent -> centroid[0]) + cell_length;
		  Vtx -> loc[1] = Parent -> centroid[1];
		  Vtx -> loc[2] = Parent -> centroid[2];
		}
	      else if(dir == NTH)
		{
		  Vtx -> loc[0] = Parent -> centroid[0];
		  Vtx -> loc[1] = (Parent -> centroid[1]) + cell_length;
		  Vtx -> loc[2] = Parent -> centroid[2];
		}
	      else
		{
		  Vtx -> loc[0] = Parent -> centroid[0];
		  Vtx -> loc[1] = Parent -> centroid[1];
		  Vtx -> loc[2] = (Parent -> centroid[2]) + cell_length;
		}
	    }
      	 
	  /*By now a vertex will exist*/

	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Assign each child on this face to this vertex*/
	    {
	      switch(dir)
		{
		case EST:
		  if(j <= 1)
		    {
		      child = j+2; /*Child on this face*/
		      mid_vertex = 7-j; /*Corresponding vertex no.*/
		    }
		  else
		    {
		      child = j+4;
		      if(j == 2)
			mid_vertex = 3;
		      else mid_vertex = 2;
		    }
		  break;
		case NTH:
		  child = 2*j+1;
		  mid_vertex = 7-2*j;
		  break;
		case UPR:
		  child = j+4;
		  mid_vertex = 7-j;
		  break;
		}

	      child_cells[child] -> verticies[mid_vertex] = Vtx;
	      Vtx -> Leaf_cells[7 - mid_vertex] = child_cells[child];
	    }
	}

      /*Now assign mid-edge verticies to children - see if neighbours already have this vertex though (this is the most
	cumbersome part)*/
      
      /*LINE SPANNING VTX 3 - VTX 7 OF PARENT*/
      
      Vtx = NULL;

      if(child_cells[7] -> face_neighbours[EST][0] != NULL) /*Try east direction*/
	{
	  if(child_cells[7] -> face_neighbours[EST][0] -> cell_level == cell_level)
	    Vtx = child_cells[7] -> face_neighbours[EST][0] -> verticies[1];
	  else if(child_cells[7] -> face_neighbours[EST][1] != NULL)
	    Vtx = child_cells[7] -> face_neighbours[EST][1] -> verticies[1];
	}
    
      if(Vtx == NULL)
	{
	  if(child_cells[7] -> face_neighbours[NTH][0] != NULL) /*Try north direction*/
	    {
	      if(child_cells[7] -> face_neighbours[NTH][0] -> cell_level == cell_level)
		Vtx = child_cells[7] -> face_neighbours[NTH][0] -> verticies[2];
	      else if(child_cells[7] -> face_neighbours[NTH][1] != NULL) 
		Vtx = child_cells[7] -> face_neighbours[NTH][1] -> verticies[2];
	    }
	}
     
      if(Vtx == NULL) /*Try 'diagonal' neighbour*/
	{ 
	  /*Derive centroid of diagnoal neighbour*/
	  diag_shift[0] = cell_length; diag_shift[1] = cell_length; diag_shift[2] = 0;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[7] -> centroid[i] + diag_shift[i];

	  /*Derive bounding box of diagnoal neighbour - cell 7 chosen as bounding box already known*/
	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[7] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[7] -> verticies[7] -> loc[i] + diag_shift[i];
	    }
	  
	  Diag_cell = find_vtx_cell(Parent -> verticies[7] -> Leaf_cells, diag_cen, diag_bbox);
	
	  if(Diag_cell != NULL) /*'Diagonal' neighbour found*/
	    Vtx = Diag_cell -> verticies[0];		

	  if(Vtx == NULL) /*Must allocate (in parallel case the vertex of the diagonal neighbour may not be allocated yet)*/
	    { 	      
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[7], 3, Vtx -> loc); /*Vertex 3 of cell 7, vertex 7 of cell 3*/
	    } 
	} 
     
      /*By now Vtx != NULL - add cells 7 and 3 to vertex array*/
      
      child_cells[7] -> verticies[3] = Vtx;
      child_cells[3] -> verticies[7] = Vtx;
      
      Vtx -> Leaf_cells[4] = child_cells[7];
      Vtx -> Leaf_cells[0] = child_cells[3];
     

      /*LINE SPANNING VTX 5 - VTX 7 OF PARENT*/

      Vtx = NULL;

      if(child_cells[7] -> face_neighbours[NTH][0] != NULL) /*Try north direction*/
	{
	  if(child_cells[7] -> face_neighbours[NTH][0] -> cell_level == cell_level) 
	    Vtx = child_cells[7] -> face_neighbours[NTH][0] -> verticies[4];
	  else if(child_cells[7] -> face_neighbours[NTH][1] != NULL) 
	    Vtx = child_cells[7] -> face_neighbours[NTH][2] -> verticies[4];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[7] -> face_neighbours[UPR][0] != NULL) /*Try upper direction*/
	    {
	      if(child_cells[7] -> face_neighbours[UPR][0] -> cell_level == cell_level) 
		Vtx = child_cells[7] -> face_neighbours[UPR][0] -> verticies[1];
	      else if(child_cells[7] -> face_neighbours[UPR][1] != NULL) 
		Vtx = child_cells[7] -> face_neighbours[UPR][1] -> verticies[1];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = 0; diag_shift[1] = cell_length; diag_shift[2] = cell_length;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[7] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[7] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[7] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[7] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[0];

	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;
	  
	      get_vtx_loc(child_cells[7], 5, Vtx -> loc); /*Vertex 5 of cell 7, vertex 7 of cell 5*/
	    }
	}

      child_cells[7] -> verticies[5] = Vtx;
      child_cells[5] -> verticies[7] = Vtx;

      Vtx -> Leaf_cells[2] = child_cells[7];
      Vtx -> Leaf_cells[0] = child_cells[5];


      /*LINE SPANNING VTX 6 - VTX 7 OF PARENT*/

      Vtx = NULL;

      if(child_cells[7] -> face_neighbours[EST][0] != NULL) /*Try east direction*/
	{
	  if(child_cells[7] -> face_neighbours[EST][0] -> cell_level == cell_level) 
	    Vtx = child_cells[7] -> face_neighbours[EST][0] -> verticies[4];
	  else if(child_cells[7] -> face_neighbours[EST][1] != NULL) 
	    Vtx = child_cells[7] -> face_neighbours[EST][2] -> verticies[4];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[7] -> face_neighbours[UPR][0] != NULL) /*Try upper direction*/
	    {
	      if(child_cells[7] -> face_neighbours[UPR][0] -> cell_level == cell_level) 
		Vtx = child_cells[7] -> face_neighbours[UPR][0] -> verticies[2];
	      else if(child_cells[7] -> face_neighbours[UPR][1] != NULL) 
		Vtx = child_cells[7] -> face_neighbours[UPR][2] -> verticies[2];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = cell_length; diag_shift[1] = 0; diag_shift[2] = cell_length;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[7] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[7] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[7] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[7] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[0];

	  if(Vtx == NULL)
	    {	      
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[7], 6, Vtx -> loc); /*Vertex 6 of cell 7, vertex 7 of cell 6*/
	    }
	}

      child_cells[7] -> verticies[6] = Vtx;
      child_cells[6] -> verticies[7] = Vtx;

      Vtx -> Leaf_cells[1] = child_cells[7];
      Vtx -> Leaf_cells[0] = child_cells[6];

    }
  
  if((phase == ALL) || (phase == NEGATIVE))
    {
      /*'NEGATIVE' PHASE*/

      /*Now assign mid-plane (of parent) verticies to children on 'negative' faces - but see if neighbours already have them*/

      for(i = 0; i < 3; i++) /*Cycle through each face*/
	{
	  dir = 2*i + 1; /*All 'negative' directions*/

	  if(Parent -> face_neighbours[dir][1] != NULL) /*4 neighbours on this face - vertex already exists*/
	    {
	      switch(dir)
		{
		case WST:
		  Vtx = child_cells[5] -> face_neighbours[dir][0] -> verticies[2]; break;
		case STH:
		  Vtx = child_cells[6] -> face_neighbours[dir][0] -> verticies[1]; break;
		case LWR:
		  Vtx = child_cells[3] -> face_neighbours[dir][0] -> verticies[4]; break;
		}
	    }
	  else /*Must allocate vertex - it doesn't exist yet*/
	    {
	      Vtx = NULL;
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(j = 0; j < MAX_NUM_CHILDREN; j++)
		Vtx -> Leaf_cells[j] = NULL; /*Initialize to NULL*/

	      /*We can get the vertex locations easily - simply shift centroid of parent to appropriate face*/

	      if(dir == WST)
		{
		  Vtx -> loc[0] = (Parent -> centroid[0]) - cell_length;
		  Vtx -> loc[1] = Parent -> centroid[1];
		  Vtx -> loc[2] = Parent -> centroid[2];
		}
	      else if(dir == STH)
		{
		  Vtx -> loc[0] = Parent -> centroid[0];
		  Vtx -> loc[1] = (Parent -> centroid[1]) - cell_length;
		  Vtx -> loc[2] = Parent -> centroid[2];
		}
	      else
		{
		  Vtx -> loc[0] = Parent -> centroid[0];
		  Vtx -> loc[1] = Parent -> centroid[1];
		  Vtx -> loc[2] = (Parent -> centroid[2]) - cell_length;
		}
	    }

	  /*By now a vertex will exist*/

	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Assign each child on this face to this vertex*/
	    {
	      switch(dir)
		{
		case WST:
		  if(j <= 1)
		    {
		      child = j;
		      mid_vertex = 5-j;
		    }
		  else
		    {
		      child = j+2;
		      if(j == 2)
			mid_vertex = 1;
		      else mid_vertex = 0;
		    }
		  break;
		case STH:
		  child = 2*j;
		  mid_vertex = 6-2*j;
		  break;
		case LWR:
		  child = j;
		  mid_vertex = 3-j;
		  break;
		}

	      child_cells[child] -> verticies[mid_vertex] = Vtx;
	      Vtx -> Leaf_cells[7 - mid_vertex] = child_cells[child];
	    }
	}


      /*LINE SPANNING VTX 0 - VTX 4 OF PARENT*/

      Vtx = NULL;

      /*Cell 0 chosen as its bounding box already known*/

      if(child_cells[0] -> face_neighbours[WST][0] != NULL) /*Try west direction*/
	{
	  if(child_cells[0] -> face_neighbours[WST][0] -> cell_level == cell_level) 
	    Vtx = child_cells[0] -> face_neighbours[WST][0] -> verticies[6];
	  else if(child_cells[0] -> face_neighbours[WST][1] != NULL) 
	    Vtx = child_cells[0] -> face_neighbours[WST][2] -> verticies[6];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[0] -> face_neighbours[STH][0] != NULL) /*Try south direction*/
	    {
	      if(child_cells[0] -> face_neighbours[STH][0] -> cell_level == cell_level) 
		Vtx = child_cells[0] -> face_neighbours[STH][0] -> verticies[5];
	      else if(child_cells[0] -> face_neighbours[STH][1] != NULL) 
		Vtx = child_cells[0] -> face_neighbours[STH][2] -> verticies[5];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = -1*cell_length; diag_shift[1] = -1*cell_length; diag_shift[2] = 0;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[0] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[0] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[0] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[0] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[7];
	  
	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[0], 4, Vtx -> loc); /*Vertex 4 of cell 0, vertex 0 of cell 4*/
	    }
	}

      child_cells[0] -> verticies[4] = Vtx;
      child_cells[4] -> verticies[0] = Vtx;

      Vtx -> Leaf_cells[3] = child_cells[0];
      Vtx -> Leaf_cells[7] = child_cells[4];


      /*LINE SPANNING VTX 0 - VTX 2 OF PARENT*/

      Vtx = NULL;

      if(child_cells[0] -> face_neighbours[STH][0] != NULL) /*Try south direction*/
	{
	  if(child_cells[0] -> face_neighbours[STH][0] -> cell_level == cell_level) 
	    Vtx = child_cells[0] -> face_neighbours[STH][0] -> verticies[3];
	  else if(child_cells[0] -> face_neighbours[STH][1] != NULL) 
	    Vtx = child_cells[0] -> face_neighbours[STH][1] -> verticies[3];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[0] -> face_neighbours[LWR][0] != NULL) /*Try lower direction*/
	    {
	      if(child_cells[0] -> face_neighbours[LWR][0] -> cell_level == cell_level) 
		Vtx = child_cells[0] -> face_neighbours[LWR][0] -> verticies[6];
	      else if(child_cells[0] -> face_neighbours[LWR][1] != NULL) 
		Vtx = child_cells[0] -> face_neighbours[LWR][2] -> verticies[6];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = 0; diag_shift[1] = -1*cell_length; diag_shift[2] = -1*cell_length;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[0] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[0] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[0] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[0] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[7];

	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[0], 2, Vtx -> loc); /*Vertex 2 of cell 0, vertex 0 of cell 2*/
	    }
	}

      child_cells[0] -> verticies[2] = Vtx;
      child_cells[2] -> verticies[0] = Vtx;

      Vtx -> Leaf_cells[5] = child_cells[0];
      Vtx -> Leaf_cells[7] = child_cells[2];


      /*LINE SPANNING VTX 0 - VTX 1 OF PARENT*/

      Vtx = NULL;

      if(child_cells[0] -> face_neighbours[WST][0] != NULL) /*Try west direction*/
	{
	  if(child_cells[0] -> face_neighbours[WST][0] -> cell_level == cell_level) 
	    Vtx = child_cells[0] -> face_neighbours[WST][0] -> verticies[3];
	  else if(child_cells[0] -> face_neighbours[WST][1] != NULL) 
	    Vtx = child_cells[0] -> face_neighbours[WST][1] -> verticies[3];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[0] -> face_neighbours[LWR][0] != NULL) /*Try lower direction*/
	    {
	      if(child_cells[0] -> face_neighbours[LWR][0] -> cell_level == cell_level) 
		Vtx = child_cells[0] -> face_neighbours[LWR][0] -> verticies[5];
	      else if(child_cells[0] -> face_neighbours[LWR][1] != NULL) 
		Vtx = child_cells[0] -> face_neighbours[LWR][1] -> verticies[5];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = -1*cell_length; diag_shift[1] = 0; diag_shift[2] = -1*cell_length;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[0] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[0] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[0] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[0] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[7];

	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);
	      
	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[0], 1, Vtx -> loc); /*Vertex 1 of cell 0, vertex 0 of cell 1*/
	    }
	}

      child_cells[0] -> verticies[1] = Vtx;
      child_cells[1] -> verticies[0] = Vtx;

      Vtx -> Leaf_cells[6] = child_cells[0];
      Vtx -> Leaf_cells[7] = child_cells[1];
      
    }

  if((phase == ALL) || (phase == POSITIVE2))
    {
      /*'SECONDARY POSITIVE' PHASE*/

      /*LINE SPANNING VTX 2 - VTX 6 OF PARENT*/

      /*Vertex 7 of cell 6 is now known*/

      Vtx = NULL;

      if(child_cells[6] -> face_neighbours[EST][0] != NULL) /*Try east direction*/
	{
	  if(child_cells[6] -> face_neighbours[EST][0] -> cell_level == cell_level) 
	    Vtx = child_cells[6] -> face_neighbours[EST][0] -> verticies[0];
	  else if(child_cells[6] -> face_neighbours[EST][1] != NULL) 
	    Vtx = child_cells[6] -> face_neighbours[EST][0] -> verticies[0];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[6] -> face_neighbours[STH][0] != NULL) /*Try south direction*/
	    {
	      if(child_cells[6] -> face_neighbours[STH][0] -> cell_level == cell_level)
		Vtx = child_cells[6] -> face_neighbours[STH][0] -> verticies[3];
	      else if(child_cells[6] -> face_neighbours[STH][1] != NULL)
		Vtx = child_cells[6] -> face_neighbours[STH][1] -> verticies[3];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = cell_length; diag_shift[1] = -1*cell_length; diag_shift[2] = 0;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[6] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[6] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[6] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[6] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[1];
	  
	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[6], 2, Vtx -> loc); /*Vertex 2 of cell 6, vertex 6 of cell 2*/
	    }
	}

      child_cells[6] -> verticies[2] = Vtx;
      child_cells[2] -> verticies[6] = Vtx;

      Vtx -> Leaf_cells[5] = child_cells[6];
      Vtx -> Leaf_cells[1] = child_cells[2];


      /*LINE SPANNING VTX 1 - VTX 3 OF PARENT*/

      /*Vertex 7 of cell 3 is known now*/

      Vtx = NULL;

      if(child_cells[3] -> face_neighbours[NTH][0] != NULL) /*Try north direction*/
	{
	  if(child_cells[3] -> face_neighbours[NTH][0] -> cell_level == cell_level) 
	    Vtx = child_cells[3] -> face_neighbours[NTH][0] -> verticies[0];
	  else if(child_cells[3] -> face_neighbours[NTH][1] != NULL) 
	    Vtx = child_cells[3] -> face_neighbours[NTH][0] -> verticies[0];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[3] -> face_neighbours[LWR][0] != NULL) /*Try lower direction*/
	    {
	      if(child_cells[3] -> face_neighbours[LWR][0] -> cell_level == cell_level) 
		Vtx = child_cells[3] -> face_neighbours[LWR][0] -> verticies[5];
	      else if(child_cells[3] -> face_neighbours[LWR][1] != NULL) 
		Vtx = child_cells[3] -> face_neighbours[LWR][1] -> verticies[5];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = 0; diag_shift[1] = cell_length; diag_shift[2] = -1*cell_length;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[3] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[3] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[3] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[3] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[4];
	  
	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[3], 1, Vtx -> loc); /*Vertex 1 of cell 3, vertex 3 of cell 1*/
	    }
	}

      child_cells[3] -> verticies[1] = Vtx;
      child_cells[1] -> verticies[3] = Vtx;

      Vtx -> Leaf_cells[6] = child_cells[3];
      Vtx -> Leaf_cells[4] = child_cells[1];


      /*LINE SPANNING VTX 2 - VTX 3 OF PARENT*/

      /*Vertex 7 of cell 3 is now known*/

      Vtx = NULL;

      if(child_cells[3] -> face_neighbours[EST][0] != NULL) /*Try east direction*/
	{
	  if(child_cells[3] -> face_neighbours[EST][0] -> cell_level == cell_level) 
	    Vtx = child_cells[3] -> face_neighbours[EST][0] -> verticies[0];
	  else if(child_cells[3] -> face_neighbours[EST][1] != NULL) 
	    Vtx = child_cells[3] -> face_neighbours[EST][0] -> verticies[0];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[3] -> face_neighbours[LWR][0] != NULL) /*Try lower direction*/
	    {
	      if(child_cells[3] -> face_neighbours[LWR][0] -> cell_level == cell_level) 
		Vtx = child_cells[3] -> face_neighbours[LWR][0] -> verticies[6];
	      else if(child_cells[3] -> face_neighbours[LWR][1] != NULL) 
		Vtx = child_cells[3] -> face_neighbours[LWR][2] -> verticies[6];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = cell_length; diag_shift[1] = 0; diag_shift[2] = -1*cell_length;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[3] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[3] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[3] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[3] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[4];
	  
	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);
	  
	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[3], 2, Vtx -> loc); /*Vertex 2 of cell 3, vertex 3 of cell 2*/
	    }
	}

      child_cells[3] -> verticies[2] = Vtx;
      child_cells[2] -> verticies[3] = Vtx;

      Vtx -> Leaf_cells[5] = child_cells[3];
      Vtx -> Leaf_cells[4] = child_cells[2];
    }

  if((phase == ALL) || (phase == NEGATIVE2))
    {
      /*'SECONDARY NEGATIVE' PHASE*/

      /*LINE SPANNING VTX 1 - VTX 5 OF PARENT*/

      /*Vertex 0 of cell 1 is now known*/

      Vtx = NULL;

      if(child_cells[1] -> face_neighbours[WST][0] != NULL) /*Try west direction*/
	{
	  if(child_cells[1] -> face_neighbours[WST][0] -> cell_level == cell_level) 
	    Vtx = child_cells[1] -> face_neighbours[WST][0] -> verticies[7];
	  else if(child_cells[1] -> face_neighbours[WST][1] != NULL) 
	    Vtx = child_cells[1] -> face_neighbours[WST][3] -> verticies[7];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[1] -> face_neighbours[NTH][0] != NULL) /*Try north direction*/
	    {
	      if(child_cells[1] -> face_neighbours[NTH][0] -> cell_level == cell_level)
		Vtx = child_cells[1] -> face_neighbours[NTH][0] -> verticies[4];
	      else if(child_cells[1] -> face_neighbours[NTH][1] != NULL) 
		Vtx = child_cells[1] -> face_neighbours[NTH][2] -> verticies[4];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = -1*cell_length; diag_shift[1] = cell_length; diag_shift[2] = 0;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[1] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[1] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[1] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[1] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[6];
	  
	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[1], 5, Vtx -> loc); /*Vertex 5 of cell 1, vertex 1 of cell 5*/
	    }
	}

      child_cells[1] -> verticies[5] = Vtx;
      child_cells[5] -> verticies[1] = Vtx;

      Vtx -> Leaf_cells[2] = child_cells[1];
      Vtx -> Leaf_cells[6] = child_cells[5];


      /*LINE SPANNING VTX 4 - VTX 6 OF PARENT*/

      Vtx = NULL;

      if(child_cells[6] -> face_neighbours[STH][0] != NULL) /*Try south direction*/
	{
	  if(child_cells[6] -> face_neighbours[STH][0] -> cell_level == cell_level) 
	    Vtx = child_cells[6] -> face_neighbours[STH][0] -> verticies[5];
	  else if(child_cells[6] -> face_neighbours[STH][1] != NULL) 
	    Vtx = child_cells[6] -> face_neighbours[STH][2] -> verticies[5];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[6] -> face_neighbours[UPR][0] != NULL) /*Try upper direction*/
	    {
	      if(child_cells[6] -> face_neighbours[UPR][0] -> cell_level == cell_level) 
		Vtx = child_cells[6] -> face_neighbours[UPR][0] -> verticies[0];
	      else if(child_cells[6] -> face_neighbours[UPR][1] != NULL) 
		Vtx = child_cells[6] -> face_neighbours[UPR][0] -> verticies[0];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = 0; diag_shift[1] = -1*cell_length; diag_shift[2] = cell_length;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[6] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[6] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[6] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[6] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[1];
	  
	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[6], 4, Vtx -> loc); /*Vertex 4 of cell 6, vertex 6 of cell 4*/
	    }
	}

      child_cells[6] -> verticies[4] = Vtx;
      child_cells[4] -> verticies[6] = Vtx;

      Vtx -> Leaf_cells[3] = child_cells[6];
      Vtx -> Leaf_cells[1] = child_cells[4];


      /*LINE SPANNING VTX 4 - VTX 5 OF PARENT*/

      /*Vertex 7 of cell 5 now known*/

      Vtx = NULL;

      if(child_cells[5] -> face_neighbours[WST][0] != NULL) /*Try west direction*/
	{
	  if(child_cells[5] -> face_neighbours[WST][0] -> cell_level == cell_level) 
	    Vtx = child_cells[5] -> face_neighbours[WST][0] -> verticies[6];
	  else if(child_cells[5] -> face_neighbours[WST][1] != NULL) 
	    Vtx = child_cells[5] -> face_neighbours[WST][2] -> verticies[6];
	}

      if(Vtx == NULL)
	{
	  if(child_cells[5] -> face_neighbours[UPR][0] != NULL) /*Try upper direction*/
	    {
	      if(child_cells[5] -> face_neighbours[UPR][0] -> cell_level == cell_level) 
		Vtx = child_cells[5] -> face_neighbours[UPR][0] -> verticies[0];
	      else if(child_cells[5] -> face_neighbours[UPR][1] != NULL) 
		Vtx = child_cells[5] -> face_neighbours[UPR][0] -> verticies[0];
	    }
	}

      if(Vtx == NULL)
	{
	  diag_shift[0] = -1*cell_length; diag_shift[1] = 0; diag_shift[2] = cell_length;
	  for(i = 0; i < 3; i++)
	    diag_cen[i] = child_cells[5] -> centroid[i] + diag_shift[i];

	  for(i = 0; i < 3; i++)
	    {
	      diag_bbox[0][i] = child_cells[5] -> verticies[0] -> loc[i] + diag_shift[i];
	      diag_bbox[1][i] = child_cells[5] -> verticies[7] -> loc[i] + diag_shift[i];
	    }

	  Diag_cell = find_vtx_cell(Parent -> verticies[5] -> Leaf_cells, diag_cen, diag_bbox);

	  if(Diag_cell != NULL) 
	    Vtx = Diag_cell -> verticies[2];
	  
	  if(Vtx == NULL)
	    {
	      Vtx = malloc(sizeof(struct vertex));
      
	      if(Vtx == NULL)
		return(ERROR);
	      else Vtx -> vtx_nums = malloc(sizeof(int)*number_of_threads);

	      if(Vtx -> vtx_nums == NULL)
		return(ERROR);

	      if(add_to_Vtxlist(Vtx, Vtxhead, Vtxtail) == ERROR)
		return(ERROR);

	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		Vtx -> Leaf_cells[i] = NULL;

	      get_vtx_loc(child_cells[5], 4, Vtx -> loc); /*Vertex 4 of cell 5, vertex 5 of cell 4*/
	    }
	}

      child_cells[5] -> verticies[4] = Vtx;
      child_cells[4] -> verticies[5] = Vtx;

      Vtx -> Leaf_cells[3] = child_cells[5];
      Vtx -> Leaf_cells[2] = child_cells[4];
    }
    
  /*All verticies have now been assigned*/

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief In serial, put parent back on vertex lists after coarsening and remove children; deallocate vertex if no longer held by any 
   cell*/

short int assign_vtx_after_coarsen(Cart_cell child_cells[], short int phase, List_vtx_glob *Vtxhead, List_vtx_glob *Vtxtail) 
{
  short int i, j, dir;
  short int pos[4]; /*Position of cells in cell array of verticies*/
  Vertex Vtx;
  Cart_cell Parent = child_cells[0] -> parent;
  Vtx = NULL;
  if((phase == ALL) || (phase == POSITIVE))
    {
      /*'POSITIVE' PHASE*/

      /*Delete volume centre vertex*/
      if(adapt_in_parallel == FALSE)
	{
	  delete_from_Vtxlist(child_cells[0] -> verticies[7], Vtxhead, Vtxtail);
	  free(child_cells[0] -> verticies[7] -> vtx_nums);
	  free(child_cells[0] -> verticies[7]); /*Children will be deallocated anyway*/
	}
      else child_cells[0] -> verticies[7] -> Vtxlist_loc -> delete = TRUE; /*If adapting in parallel, just flag for deletion*/

      /*Add parent back to corner verticies*/
      for(i = 0; i < MAX_NUM_CHILDREN; i++) 
	Parent -> verticies[i] -> Leaf_cells[7-i] = Parent;

      /*Delete face verticies if necessary*/
      for(i = 0; i < 3; i++)
	{
	  dir = 2*i;

	  switch(dir)
	    {
	    case EST:
	      Vtx = child_cells[2] -> verticies[7]; /*Planar vertex on this parental face*/
	      pos[0] = 0; /*Cells on east parental face occupy these positions on the centre vertex there*/ 
	      pos[1] = 1; 
	      pos[2] = 4;
	      pos[3] = 5;
	      break;
	    case NTH:
	      Vtx = child_cells[1] -> verticies[7];
	      pos[0] = 0;
	      pos[1] = 2;
	      pos[2] = 4;
	      pos[3] = 6;
	      break;
	    case UPR:
	      Vtx = child_cells[4] -> verticies[7];
	      pos[0] = 0;
	      pos[1] = 1;
	      pos[2] = 2;
	      pos[3] = 3;
	      break;
	    }

	  for(j = 0; j < 4; j++)
	    Vtx -> Leaf_cells[pos[j]] = NULL;

	  for(j = 0; j < MAX_NUM_CHILDREN; j++)
	    {
	      if(Vtx -> Leaf_cells[j] != NULL)
		break;
	    }

	  if(j == MAX_NUM_CHILDREN) /*No cell needs this vertex anymore*/
	    {
	      if(adapt_in_parallel == FALSE)
		{
		  delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
		  free(Vtx -> vtx_nums);
		  free(Vtx);
		}
	      else Vtx -> Vtxlist_loc -> delete = TRUE;
	    }
	}

      /*Now free edge verticies if necessary*/

      /*VERTEX BETWEEN VERTICIES 3 AND 7*/
      Vtx = child_cells[3] -> verticies[7];
      Vtx -> Leaf_cells[0] = NULL; 
      Vtx -> Leaf_cells[4] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}

      /*VERTEX BETWEEN VERTICIES 5 AND 7*/
      Vtx = child_cells[5] -> verticies[7];
      Vtx -> Leaf_cells[0] = NULL; 
      Vtx -> Leaf_cells[2] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}

      /*VERTEX BETWEEN VERTICIES 6 AND 7*/
      Vtx = child_cells[6] -> verticies[7];
      Vtx -> Leaf_cells[0] = NULL; 
      Vtx -> Leaf_cells[1] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}
    }
  
  /*'NEGATIVE' PHASE*/

  if((phase == ALL) || (phase == NEGATIVE))
    {
      /*Delete face verticies if necessary*/
      for(i = 0; i < 3; i++)
	{
	  dir = 2*i + 1;

	  switch(dir)
	    {
	    case WST:
	      Vtx = child_cells[5] -> verticies[0]; 
	      pos[0] = 2;
	      pos[1] = 3;
	      pos[2] = 6;
	      pos[3] = 7;
	      break;
	    case STH:
	      Vtx = child_cells[0] -> verticies[6];
	      pos[0] = 1;
	      pos[1] = 3;
	      pos[2] = 5;
	      pos[3] = 7;
	      break;
	    case LWR:
	      Vtx = child_cells[3] -> verticies[0];
	      pos[0] = 4;
	      pos[1] = 5;
	      pos[2] = 6;
	      pos[3] = 7; 
	      break;
	    }

	  for(j = 0; j < 4; j++)
	    Vtx -> Leaf_cells[pos[j]] = NULL;

	  for(j = 0; j < MAX_NUM_CHILDREN; j++)
	    {
	      if(Vtx -> Leaf_cells[j] != NULL)
		break;
	    }

	  if(j == MAX_NUM_CHILDREN) /*No cell needs this vertex anymore*/
	    {
	      if(adapt_in_parallel == FALSE)
		{
		  delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
		  free(Vtx -> vtx_nums);
		  free(Vtx);
		}
	      else Vtx -> Vtxlist_loc -> delete = TRUE;
	    }
	}

      /*VERTEX BETWEEN VERTICIES 0 AND 4*/
      Vtx = child_cells[0] -> verticies[4];
      Vtx -> Leaf_cells[3] = NULL; 
      Vtx -> Leaf_cells[7] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}

      /*VERTEX BETWEEN VERTICIES 0 AND 2*/
      Vtx = child_cells[0] -> verticies[2];
      Vtx -> Leaf_cells[5] = NULL; 
      Vtx -> Leaf_cells[7] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}

      /*VERTEX BETWEEN VERTICIES 0 AND 1*/
      Vtx = child_cells[0] -> verticies[1];
      Vtx -> Leaf_cells[6] = NULL; 
      Vtx -> Leaf_cells[7] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}
    } 

  /*'SECONDARY POSITIVE' PHASE*/

  if((phase == ALL) || (phase == POSITIVE2))
    {
      /*VERTEX BETWEEN VERTICIES 2 AND 6*/
      Vtx = child_cells[2] -> verticies[6];
      Vtx -> Leaf_cells[1] = NULL; 
      Vtx -> Leaf_cells[5] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}

      /*VERTEX BETWEEN VERTICIES 1 AND 3*/
      Vtx = child_cells[1] -> verticies[3];
      Vtx -> Leaf_cells[4] = NULL; 
      Vtx -> Leaf_cells[6] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}

      /*VERTEX BETWEEN VERTICIES 2 AND 3*/
      Vtx = child_cells[2] -> verticies[3];
      Vtx -> Leaf_cells[4] = NULL; 
      Vtx -> Leaf_cells[5] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}
    }

  /*'SECONDARY NEGATIVE' PHASE*/

  if((phase == ALL) || (phase == NEGATIVE2))
    {
      /*VERTEX BETWEEN VERTICIES 1 AND 5*/
      Vtx = child_cells[1] -> verticies[5];
      Vtx -> Leaf_cells[2] = NULL; 
      Vtx -> Leaf_cells[6] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}

      /*VERTEX BETWEEN VERTICIES 4 AND 6*/
      Vtx = child_cells[4] -> verticies[6];
      Vtx -> Leaf_cells[1] = NULL; 
      Vtx -> Leaf_cells[3] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}

      /*VERTEX BETWEEN VERTICIES 4 AND 5*/
      Vtx = child_cells[4] -> verticies[5];
      Vtx -> Leaf_cells[2] = NULL; 
      Vtx -> Leaf_cells[3] = NULL;
      for(j = 0; j < MAX_NUM_CHILDREN; j++)
	{
	  if(Vtx -> Leaf_cells[j] != NULL)
	    break;
	}
      if(j == MAX_NUM_CHILDREN)
	{
	  if(adapt_in_parallel == FALSE)
	    {
	      delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	      free(Vtx -> vtx_nums);
	      free(Vtx);
	    }
	  else Vtx -> Vtxlist_loc -> delete = TRUE;
	}
    }
  
  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Count the total no. of verticies and split them between threads*/

void count_verticies(void)
{
  int num_vtxs = 0;
  List_vtx_glob L;
  int vtxs_per_thread; 
  int vtxs = 1; /*Counter for no. cells in a sublist*/
  int thread_num = 0; 
  
  L = Vtxlist;
  while(L != NULL)
    {
      num_vtxs++;
      L = L -> next;
    }

  /*Get no. of threads for each cell (except last thread)*/
  vtxs_per_thread = ROUND(((double) num_vtxs)/((double) number_of_threads));

  if(num_vtxs <= number_of_threads)
    {
      /*Unlikely?  But reduce no. threads so num_vtxs == number_of_threads, then must allocate new sublist pointers*/
    }

  /*Now traverse list and vtxs as belonging to different sublists*/
  
  L = Vtxlist; 
  Vtx_heads[0] = L; 
  while(L != NULL)
    {
      if((vtxs <= vtxs_per_thread) || (thread_num == (number_of_threads - 1)))
	{	
	  if((vtxs == vtxs_per_thread) || (L -> next == NULL)) 
	    Vtx_tails[thread_num] = L; /*Tail of sublist*/	  
	  vtxs++;  
	}
      else 
	{
	  vtxs = 1; /*Restart counter*/
	  thread_num++; /*Increase thread number*/
	  Vtx_heads[thread_num] = L; /*Head of sublist*/
	  vtxs++;
	}

      L = L -> next;
    }
}

/*------------------------------------------------------------------*/

/**\brief Helper function for vertex assigning - if a cell on the cell array has the required
   centroid or its descendants contained within its bounding box, will return the cell with the centroid*/

Cart_cell find_vtx_cell(Cart_cell Leaf_cells[], double centroid[], double bbox[][3])
{
  short int i, j;
  Cart_cell C;
    
  for(i = 0; i < MAX_NUM_CHILDREN; i++) /*Go through array of leaf cells*/
    {
      C = Leaf_cells[i];

      if(C != NULL)
	{
	  if(
#if FLOAT_EQ_STABLE
	     EQ(C -> centroid[0], centroid[0]) && EQ(C -> centroid[1], centroid[1]) && EQ(C -> centroid[2], centroid[2])
	     /*Shouldn't have cells too small, else these comparisons are stuffed*/
#else
	     (C -> centroid[0] == centroid[0]) && (C -> centroid[1] == centroid[1]) && (C -> centroid[2] == centroid[2])
#endif
	     )
	    return(C);
	  else
	    { /*Check if cell is inside bounding box - if so one of its ancestors will have the right centroid*/
	      
	      for(j = 0; j < 3; j++)
		{
		  if((((C -> centroid[j]) - 0.5*CALC_CELL_EDGE_LENGTH(C)) < bbox[0][j]) || 
		     (((C -> centroid[j]) + 0.5*CALC_CELL_EDGE_LENGTH(C)) > bbox[1][j]))
		    break;
		}  
	      
	      if(j == 3) /*Cell is within bounding box*/
		{ /*Keep traversing the ancestors until the correct cell is found*/

		  while((C != NULL) && 
#if FLOAT_EQ_STABLE
			(NOT_EQ(C->centroid[0], centroid[0]) || NOT_EQ(C->centroid[1], centroid[1]) || NOT_EQ(C->centroid[2], centroid[2]))
#else
			((C->centroid[0] != centroid[0]) || (C->centroid[1] != centroid[1]) || (C->centroid[2] != centroid[2]))
#endif
			)
		    C = C -> parent;
                                                                                                                            
		  return(C);
		}
	    }
	}
    }
  
  return(NULL); /*No cell can be found*/
}

/*------------------------------------------------------------------*/

/**\brief Allocate flux vectors for cell interfaces (only need to do this for adapted cells before each timestep)*/

short int alloc_fluxes(List_leaf L, List_leaf Tail, short int directions)
{ 
  /*Because flux vector allocation only performed on adapted cells, these cells wouldn't be merged in any direction*/

  List_leaf node = L;
  Cart_cell C;
  short int i, j, dir, opp_dir;
  short int dirs[3];
  short int num_faces;
  j = 0; opp_dir = 0; num_faces = 0;
  if(directions == ALL)
    num_faces = NUM_FACES;
  else if(directions == POSITIVE)
    {
      dirs[0] = EST; dirs[1] = NTH; dirs[2] = UPR;
      num_faces = 3;
    }
  else if(directions == NEGATIVE)
    {
      dirs[0] = WST; dirs[1] = STH; dirs[2] = LWR;
      num_faces = 3;
    }

  while(node != NULL)
    {
      C = node -> cell_loc;
      
      if((C -> just_created == TRUE)  || ((C -> is_small_cell == TRUE) && (C -> Merge_data == NULL))) 
 	{ /*Only need to worry about new cells on the list*/

#if EXHAUSTIVE_DEBUG

#else
	  if((directions == ALL) || (directions == NEGATIVE))
	    C -> just_created = FALSE; /*Reset*/
#endif

	  for(dir = 0; dir < num_faces; dir++) /*Go through faces and see if any flux vectors need to be allocated*/
	    { 
	      if(num_faces == 3)
		i = dirs[dir];
	      else i = dir;

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

	      if(((C -> face_neighbours[i][0] == NULL) || (C -> face_neighbours[i][1] == NULL)) &&
		 ((C -> flux_area[i][0] > 0) || (C -> flux_area[i][1] > 0) ||
		  (C -> flux_area[i][2] > 0) || (C -> flux_area[i][3] > 0)))
		{ /*At border or only 1 neighbour and can flux through face - only 1 flux vector needed*/
#if !NOTINCLUDE
		  if((C -> Flow_data.Face_fluxes[i][0] == NULL) && 
		     ((C -> Merge_data == NULL) || (C -> Merge_data != C -> face_neighbours[i][0] -> Merge_data))) /*No need to test for
														     merged cells?*/
#else
		    if(C -> Flow_data.Face_fluxes[i][0] == NULL)
#endif 
		      { /*No flux vector allocated (it might already have been), can allocate flux vector*/
 
			C -> Flow_data.Face_fluxes[i][0] = malloc(sizeof(struct flux_vector));

			if(C -> Flow_data.Face_fluxes[i][0] == NULL)
			  return(ERROR);

			if(C -> face_neighbours[i][0] != NULL) 
			  { /*Neighbour exists - it can be same level or 1 level coarser*/

			    if((C -> face_neighbours[i][0] -> cell_level) == (C->cell_level)) /*Same level neighbour - point to same flux*/
			      C -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][0] = C -> Flow_data.Face_fluxes[i][0];
			    else /*1 level coarser*/
			      { /*Must know corresponding quadrant of neighbour*/
			  
				switch(i)
				  {
				  case NTH:
				    j = ((C -> child_num)-1)/2; break;
				    /*jth quadrant of parent cell that child resides on; neighbour's quadrant also jth quadrant*/
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

				C -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][j] = C -> Flow_data.Face_fluxes[i][0];
			      }
			  }		      
		      } /*Flux vector already exists*/
		}
	      else if(C -> face_neighbours[i][1] != NULL) /*4 cells in this direction; may have to allocate/use more than 1 flux vector*/
		{
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Cycle through each quadrant*/
		    {
		      if(C -> flux_area[i][j] > 0) /*Can flux through quadrant*/
			{  
#if !NOTINCLUDE
			  if((C -> Flow_data.Face_fluxes[i][j] == NULL) && 
			     ((C -> Merge_data == NULL) || (C -> face_neighbours[i][j] -> Merge_data != C -> Merge_data))) 
#else
			    if(C -> Flow_data.Face_fluxes[i][j] == NULL)
#endif
			      { 			      
				C -> Flow_data.Face_fluxes[i][j] = malloc(sizeof(struct flux_vector));
			  
				if(C -> Flow_data.Face_fluxes[i][j] == NULL)
				  return(ERROR);
			      			  
				C -> face_neighbours[i][j] -> Flow_data.Face_fluxes[opp_dir][0] = C -> Flow_data.Face_fluxes[i][j];
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

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief After parent is refined, childrens' faces must now point to corresponding interface
   flux vector that parent previously pointed to.*/

void reassign_fluxes_after_refine(Cart_cell child_cells[], short int directions)
{
  short int i, j, child, dir, opp_dir;
  Cart_cell Parent = child_cells[0] -> parent;

  short int dirs[3];
  short int num_faces;
  child = 0; opp_dir = 0; num_faces = 0;
  if(directions == ALL)
    num_faces = NUM_FACES;
  else if(directions == POSITIVE)
    {
      dirs[0] = EST; dirs[1] = NTH; dirs[2] = UPR;
      num_faces = 3;
    }
  else if(directions == NEGATIVE)
    {
      dirs[0] = WST; dirs[1] = STH; dirs[2] = LWR;
      num_faces = 3;
    }    

  if(directions != NEGATIVE)
    {
      if(Parent -> Flow_data.Obstructed_flux != NULL) /*Find first child that has face obstruction and reassign it*/
	{
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((child_cells[i] -> cell_type != SOLID) && (child_cells[i] -> cell_type != UNKNOWN)) /*Only give to known leaf cells*/
		{
		  for(j = 0; j < 3; j++) /*Check if net obstruction along axis (net difference in opposite areas)*/
		    {
		      dir = 2*j; /*'Positive' directions*/
		      opp_dir = dir+1; /*'Negative' directions*/

		      if(ABS((child_cells[i]->flux_area[opp_dir][0]+child_cells[i]->flux_area[opp_dir][1]+
			      child_cells[i]->flux_area[opp_dir][2]+child_cells[i]->flux_area[opp_dir][3]) - 
			     (child_cells[i]->flux_area[dir][0]+child_cells[i]->flux_area[dir][1]+
			      child_cells[i]->flux_area[dir][2]+child_cells[i]->flux_area[dir][3])) > 0)
			break;
		    }

		  if(j < 3)
		    {
		      child_cells[i] -> Flow_data.Obstructed_flux = Parent -> Flow_data.Obstructed_flux; /*Vector reassigned*/
		      break;
		    }
		}
	    }	
	}
    }
    
  /*Now inter-cell fluxes - need to get corresponding child on ith face, jth quadrant of parent*/

  for(dir = 0; dir < num_faces; dir++) /*Cycle through each face/direction*/
    {      
      if(num_faces == 3)
	i = dirs[dir];
      else i = dir;

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

      if(((Parent -> face_neighbours[i][0] == NULL) || ((Parent -> face_neighbours[i][1] == NULL))) && 
	 ((Parent -> flux_area[i][0] > 0) || (Parent -> flux_area[i][1] > 0) ||
	  (Parent -> flux_area[i][2] > 0) || (Parent -> flux_area[i][3] > 0))) /*1 flux vector in this direction stored*/
	{ /*All children have this neighbour too (which will be same level as parent) - find the first child
	    which needs this flux vector and reassign it*/
	 
	  if((Parent -> face_neighbours[i][0] != NULL) && (Parent -> face_neighbours[i][0] -> cell_level == (Parent -> cell_level-1)))
	    { /*Parent can flux in this direction but sees a coarser neighbour - just delete the flux vector as this
		case occurs when refining in parallel*/

	      if(Parent -> Flow_data.Face_fluxes[i][0] != NULL)
		{
		  free(Parent -> Flow_data.Face_fluxes[i][0]);
		  Parent -> Flow_data.Face_fluxes[i][0] = NULL;
		}
	    }
	  else
	    {
   	      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		{
		  switch(i)
		    {
		    case NTH:
		      child = 2*j+1; break;
		    case STH:
		      child = 2*j; break;
		    case EST:
		      if(j <= 1)
			child = j+2;
		      else child = j+4;
		      break;
		    case WST:
		      if(j <= 1)
			child = j;
		      else child = j+2;
		      break;
		    case LWR:
		      child = j; break;
		    case UPR:
		      child = j+4; break;
		    }

		  if((child_cells[child] -> flux_area[i][0] > 0) || (child_cells[child] -> flux_area[i][1] > 0) ||
		     (child_cells[child] -> flux_area[i][2] > 0) || (child_cells[child] -> flux_area[i][3] > 0)) /*Flux in this direction*/
		    {		  
		      child_cells[child] -> Flow_data.Face_fluxes[i][0] = (Parent -> Flow_data.Face_fluxes[i][0]);
		      Parent -> Flow_data.Face_fluxes[i][0] = NULL; /*Make it NULL or it might be later left 'hanging' with nothing to 
								      point to*/

		      if(Parent -> face_neighbours[i][0] != NULL) /*Let flux vector be the 'jth' flux vector (in direction i) for 
								    the neighbour*/
			{
			  Parent -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][0] = NULL; /*Previously it was the 0th flux
													 vector for the neighbour*/
			  Parent -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][j] = 
			    (child_cells[child] -> Flow_data.Face_fluxes[i][0]);
			}

		      break;		  
		    }
		}

	      if(j == MAX_NUM_NEIGHBOURS) /*Parent can flux thru face but no child can because of volume degeneracy - delete flux vector*/
		{
		  if(Parent -> face_neighbours[i][0] != NULL) 
		    Parent -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][0] = NULL;

		  if(Parent -> Flow_data.Face_fluxes[i][0] != NULL)
		    free(Parent -> Flow_data.Face_fluxes[i][0]);

		  Parent -> Flow_data.Face_fluxes[i][0] = NULL;
		}	    
	    }
	}
      else if(Parent -> face_neighbours[i][1] != NULL)
	{
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
	    {
	      switch(i)
		{
		case NTH:
		  child = 2*j+1; break;
		case STH:
		  child = 2*j; break;
		case EST:
		  if(j <= 1)
		    child = j+2;
		  else child = j+4;
		  break;
		case WST:
		  if(j <= 1)
		    child = j;
		  else child = j+2;
		  break;
		case LWR:
		  child = j; break;
		case UPR:
		  child = j+4; break;
		}

	      if((Parent -> flux_area[i][j] > 0) && (Parent -> face_neighbours[i][j] -> children[0] != NULL))
		Parent -> Flow_data.Face_fluxes[i][j] = NULL; /*Parent can flux thru this quadrant, but its neighbour has children - 
								ignore this case (as it occurs in parallel)*/
	      else
		{
		  if((child_cells[child] -> flux_area[i][0] > 0) || (child_cells[child] -> flux_area[i][1] > 0) ||
		     (child_cells[child] -> flux_area[i][2] > 0) || (child_cells[child] -> flux_area[i][3] > 0))  
		    {		  
		      child_cells[child] -> Flow_data.Face_fluxes[i][0] = (Parent -> Flow_data.Face_fluxes[i][j]);
		      Parent -> Flow_data.Face_fluxes[i][j] = NULL;
		    } 
		  else if(Parent -> Flow_data.Face_fluxes[i][j] != NULL) /*Again volume degenerate case (parent can flux through
									   this quad but children can't)*/
		    {
		      Parent -> face_neighbours[i][j] -> Flow_data.Face_fluxes[opp_dir][0] = NULL;
		      free(Parent -> Flow_data.Face_fluxes[i][j]);
		      Parent -> Flow_data.Face_fluxes[i][j] = NULL;
		    }
		}	      
	    }
	}
    }
}

/*------------------------------------------------------------------*/

/**\brief After parent is coarsened, it can use its childrens' interface flux vectors (by reassigning), or might
   need to deallocate them if some interfaces no longer exist.  This is essentially the inverse of the 
   reassign_fluxes_after_refine() operation*/

void reassign_fluxes_after_coarsen(Cart_cell child_cells[], short int directions)
{
  short int i, j, child, dir, opp_dir, child_found;
  short int Obstructed_flux_set = UNKNOWN;
  Cart_cell Parent = child_cells[0] -> parent;

  double parent_quad_area = CALC_CELL_FLUID_AREA(Parent)/MAX_NUM_NEIGHBOURS;

  short int dirs[3] = {EST, NTH, UPR};
  short int num_faces;
  child = 0; opp_dir = 0; num_faces = 0;
  if(directions == ALL)
    num_faces = NUM_FACES;
  else if(directions == POSITIVE)
    num_faces = 3;
  else if(directions == NEGATIVE)
    {
      dirs[0] = WST; dirs[1] = STH; dirs[2] = LWR;
      num_faces = 3;
    }    

  if(directions != NEGATIVE)
    {
      /*First reassign Obstructed_flux vector if necessary*/

      if(Parent -> Flow_data.Obstructed_flux != NULL) /*Parent may already have this flux vector*/
	Obstructed_flux_set = TRUE;
      else
	{
	  for(i = 0; i < NUM_FACES; i++)
	    {
	      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		{
		  if(Parent -> flux_area[i][j] < parent_quad_area) /*Obstruction in this direction*/
		    {
		      Obstructed_flux_set = FALSE; /*Need to find child*/
		      break;
		    } 
		}
	  
	      if(Obstructed_flux_set == FALSE)
		break;
	    }
	}

      if(Obstructed_flux_set != UNKNOWN)
	{            
	  for(i = 0; i < MAX_NUM_CHILDREN; i++) /*Test if child has this vector - either reassign or deallocate*/
	    {
	      if(child_cells[i] -> Flow_data.Obstructed_flux != NULL)
		{
		  if(Obstructed_flux_set == FALSE) /*Use a child's vector if parent not already have*/
		    {
		      Parent -> Flow_data.Obstructed_flux = child_cells[i] -> Flow_data.Obstructed_flux;
		      Obstructed_flux_set = TRUE;
		    }
		  else /*Deallocate other Obstructed_flux vectors*/
		    {
		      if(child_cells[i] -> Flow_data.Obstructed_flux != Parent -> Flow_data.Obstructed_flux)
			free(child_cells[i] -> Flow_data.Obstructed_flux);
		  		  
		      /*Parent may already point to this flux vector - take care not to deallocate it*/
		    }
		}
	    }
	}

      /*Remove any inter-sibling flux vectors first (need only cycle through positive internal faces)*/

      for(dir = 0; dir < 3; dir++)
	{	  
	  i = dirs[dir];
	 
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
	    {
	      switch(i)
		{
		case EST:
		  if(j <= 1)
		    child = j;
		  else child = j+2;
		  break;
		case NTH:
		  child = 2*j; break;
		case UPR:
		  child = j; break;
		}

	      if(child_cells[child] -> Flow_data.Face_fluxes[i][0] != NULL)
		free(child_cells[child] -> Flow_data.Face_fluxes[i][0]);
	    } 
	}
    }

  /*Now do inter-cell fluxes*/

  for(dir = 0; dir < num_faces; dir++) /*Cycle through each face of parent*/
    {
      if(num_faces == 3)
	i = dirs[dir];
      else i = dir;

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

      if(((Parent -> face_neighbours[i][0] == NULL) || (Parent -> face_neighbours[i][1] == NULL)) &&
	 ((Parent -> flux_area[i][0] > 0) || (Parent -> flux_area[i][1] > 0) ||
	  (Parent -> flux_area[i][2] > 0) || (Parent -> flux_area[i][3] > 0)))
									
	{  /*1 flux vector needed, deallocate others - if it can't flux through the face, no vectors would exist already.  Find
	     the first child that has the flux vector and use it*/
	     
	  child_found = FALSE;

	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++) /*Identify child cells on this face*/
	    {
	      switch(i)
		{
		case NTH:
		  child = 2*j+1; break;
		case STH:
		  child = 2*j; break;
		case EST:
		  if(j <= 1)
		    child = j+2;
		  else child = j+4;
		  break;
		case WST:
		  if(j <= 1)
		    child = j;
		  else child = j+2;
		  break;
		case LWR:
		  child = j; break;
		case UPR:
		  child = j+4; break;
		}

	      if((Parent -> face_neighbours[i][0] != NULL) && (Parent -> face_neighbours[i][0] -> cell_level == (Parent -> cell_level-1)))
		{
		  if(child_cells[child] -> Flow_data.Face_fluxes[i][0] != NULL)
		    free(child_cells[child] -> Flow_data.Face_fluxes[i][0]); 

		  /*Parent can flux in this direction yet neighbour is coarser - free flux vectors (this occurs in parallel adaptation)*/
		}
	      else if(child_cells[child] -> Flow_data.Face_fluxes[i][0] != NULL)
		{ 
		  if(child_found == FALSE) /*Reassign this vector to the parent*/ 
		    {		  
		      Parent -> Flow_data.Face_fluxes[i][0] = child_cells[child] -> Flow_data.Face_fluxes[i][0];
		      
		      /*No need to make child's pointer to this vector NULL; if cell refined it will be reassigned back to it*/
		      
		      if(Parent -> face_neighbours[i][0] != NULL) /*Same level neighbour*/
			{
			  Parent -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][j] = NULL; /*Previously it was the jth
													 vector for the neighbour*/

			  Parent -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][0] = 
			    Parent -> Flow_data.Face_fluxes[i][0];			
			} 
		      
		      child_found = TRUE;		  
		    }
		  else /*Can deallocate - no need to set pointers of children to NULL as they're deallocated anyway*/
		    { 
		      free(child_cells[child] -> Flow_data.Face_fluxes[i][0]);
		      /*child_cells[child] -> Flow_data.Face_fluxes[i][0] = NULL;*/
		      
		      if(Parent -> face_neighbours[i][0] != NULL)
			Parent -> face_neighbours[i][0] -> Flow_data.Face_fluxes[opp_dir][j] = NULL;
		    }
		}	          
	    }
	    	  
	}
      else if(Parent -> face_neighbours[i][1] != NULL)
	{ /*Childrens' flux vector pointers can simply be reassigned to quadrants of parent even if they're NULL*/
	  
	  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
	    {
	      switch(i)
		{
		case NTH:
		  child = 2*j+1; break;
		case STH:
		  child = 2*j; break;
		case EST:
		  if(j <= 1)
		    child = j+2;
		  else child = j+4;
		  break;
		case WST:
		  if(j <= 1)
		    child = j;
		  else child = j+2;
		  break;
		case LWR:
		  child = j; break;
		case UPR:
		  child = j+4; break;
		}
	      
	      if((Parent -> flux_area[i][j] > 0) && 
		 (child_cells[child] -> face_neighbours[i][0] -> cell_level == (child_cells[child] -> cell_level + 1)))
		{ /*Parent can flux through this quadrant and yet it has more than 1 neighbour there - don't do anything*/
		}
	      else Parent -> Flow_data.Face_fluxes[i][j] = child_cells[child] -> Flow_data.Face_fluxes[i][0];
	    }
	}
    }   
}

/*------------------------------------------------------------------*/
