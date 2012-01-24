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

/**\file 
   Source file on routines for octree cartesian VCE cells which call geometry interrogation functions - 
   computation of basis carteisan cell spatial properties, partial area/volume determinations and additional_refine_times*/

#include <stdio.h>
#include <math.h>
#include "ov_kernel.h"
#include "ov_setgeom.h"
#include "ov_adts.h"

#define NOTINCLUDE 0
#define FINE_IC 0

extern short int twoD;
extern double twoD_thickness;
extern short int use_geom_engine_always;
extern short int level_before_can_interrogate;
extern short int refining_root_now;
extern short int allow_multiple_refinements;
extern short int area_subcell_num;
extern short int volume_subcell_num;
extern short int wall_rep;

#if NOTINCLUDE
extern double IC_bbox[2][3];
#else
extern Body *IC_regions;
extern double ***IC_bboxes;
extern int num_ICs;
extern Deton_pts Deton;
#endif

extern Point_tnode ADT2Dpoint;
extern Point_tnode ADT3Dpoint;
extern short int num_times_refine_root;
extern short int min_refinement_level;
extern short int min_intersect_refine_level;
extern short int allow_multiple_refinements;
extern double root_shift[3];
extern double root_scale;

extern Body_tnode ADTbody;

/*------------------------------------------------------------------*/

/**\brief Compute centroids of cell's children relative to global co-ord sys*/

void compute_children_centroids(Cart_cell C, double centroids[][3]) 
{                                                    
  short int i, j;
  short int shift[3]; /*Translation factors to obtain children centroids*/  

  for(i = 0; i < MAX_NUM_CHILDREN; i++) /*Cycle through each child - translation factors different*/
    {
      switch(i) 
	{
	case 0: /*i.e. child 0*/
	  shift[0] = 0; shift[1] = 0; shift[2] = 0;
	  break;
	case 1:
	  shift[0] = 0; shift[1] = 1; shift[2] = 0;
	  break;
	case 2:
	  shift[0] = 1; shift[1] = 0; shift[2] = 0;
	  break;
	case 3:
	  shift[0] = 1; shift[1] = 1; shift[2] = 0;
	  break;
	case 4:
	  shift[0] = 0; shift[1] = 0; shift[2] = 1;
	  break;
	case 5:
	  shift[0] = 0; shift[1] = 1; shift[2] = 1;
	  break;
	case 6:
	  shift[0] = 1; shift[1] = 0; shift[2] = 1;
	  break;
	case 7:
	  shift[0] = 1; shift[1] = 1; shift[2] = 1;
	  break;
	}      

      for(j = 0; j < 3; j++)
	centroids[i][j] = (C -> centroid[j]) + (2*(shift[j])-1)*0.25*(C -> cell_length);
      
      /*Centroids of all children computed*/
    }
}

/*------------------------------------------------------------------*/

/**\brief Get location of vertex given by vertex no. vtx_num*/

void get_vtx_loc(Cart_cell C, short int vtx_num, double loc[]) 
{  
  short int Xt, Yt, Zt; /*Translation terms (shift from centroid to form verticies)*/
  Xt = 0; Yt = 0; Zt = 0;
  switch(vtx_num)
    {
    case 0:
      Xt = 0; Yt = 0; Zt = 0;
      break;
    case 1:
      Xt = 0; Yt = 1; Zt = 0;
      break;
    case 2:
      Xt = 1; Yt = 0; Zt = 0;
      break;
    case 3:
      Xt = 1; Yt = 1; Zt = 0;
      break;
    case 4:
      Xt = 0; Yt = 0; Zt = 1;
      break;
    case 5:
      Xt = 0; Yt = 1; Zt = 1;
      break;
    case 6:
      Xt = 1; Yt = 0; Zt = 1;
      break;
    case 7:
      Xt = 1; Yt = 1; Zt = 1;
      break;
    }

  loc[0] = (C -> centroid[0]) + (2.0*Xt-1.0)*((C -> cell_length)/2.0);
  loc[1] = (C -> centroid[1]) + (2.0*Yt-1.0)*((C -> cell_length)/2.0);
  loc[2] = (C -> centroid[2]) + (2.0*Zt-1.0)*((C -> cell_length)/2.0);

  /*Location of cell centroid + translation terms give vertex locations*/
}

/*------------------------------------------------------------------*/

/**\brief Routine to compute cell types/volumes, partial face areas and intersection coarsenabilities for
   8 children cells*/

short int set_cell_geom_data(Cart_cell child_cells[], short int refine_steps, short int do_volume, short int do_coarsenability,
			     short int do_area, short int directions, short int rm_temp_lists)
{ 
  /*refine_steps is no. times parent of child cells is refined, additionally_refine_for_geometry specifies if parent of child cells
    wants to be further refined for intersection reasons*/

  short int child_level = child_cells[0] -> cell_level; /*Level for children is the same*/
  short int i, j, k;
  double cl;
  
  if((child_level < level_before_can_interrogate) && (refining_root_now == TRUE)) {
    
    for(i = 0; i < MAX_NUM_CHILDREN; i++)
      child_cells[i] -> additional_refine_times = 1;

    if(twoD != 0) { /*Don't refine any more than we need to*/
      cl = CALC_CELL_EDGE_LENGTH(child_cells[0]);

      for(i = 0; i < 8; i++) {
	if(
#if FLOAT_EQ_STABLE
	   NOT_EQ(((child_cells[i]->centroid[2])-(cl/2.0)), (root_shift[2] - (root_scale/2.0)))
#else
	   ((child_cells[i] -> centroid[2]) - (cl/2.0)) != (root_shift[2] - (root_scale/2.0))
#endif
	   ) {
	  child_cells[i] -> cell_type = SOLID;
#if 0
	  child_cells[i] -> cell_volume = 0;
	  child_cells[i]->additional_refine_times = 0;
	  for(j = 0; j < NUM_FACES; j++)
	    for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
	      child_cells[i] -> flux_area[j][k] = 0;
#endif
	}
      }
    }
    
    return(NOERROR); /*Interrogate geometry engine only when cell is at certain level of refinement - only when root is refined,
		       else try and set geometric data anyway*/
  }
  else
    {
      double child_edge_length = CALC_CELL_EDGE_LENGTH(child_cells[0]);
      double child_face_area;
      short int face_set; /*If TRUE, face's flux area has been set*/
      short int nsdirection[3]; /*Non-sibling directions*/
      short int quadrant_num; /*The quadrant of the parent corresponding to child cell*/
      short int opp_dir; /*Opposite direction*/

      short int dir;
      short int dirs[3];
      short int num_faces;

      Cart_cell Parent = child_cells[0] -> parent; /*Parent for children is the same*/
      short int Parent_type = Parent -> cell_type; 
     
      quadrant_num = 0; opp_dir = 0; num_faces = 0; 
      /*\b Proceed to set the cell type of the children (solid, fluid or intersected) and their fluid volumes*/          
      
	if(do_volume == TRUE)
	  {
	    double child_fluid_volume = CALC_CELL_FLUID_VOL(child_cells[0]);

	    if(twoD != 0) { /*Interrogating volumes are irrelevant in 2D case*/

	      for(i = 0; i < MAX_NUM_CHILDREN; i++) {
		if(
#if FLOAT_EQ_STABLE
		   NOT_EQ(((child_cells[i]->centroid[2])-(child_edge_length/2.0)), (root_shift[2] - (root_scale/2.0)))
#else
		   ((child_cells[i] -> centroid[2]) - (child_edge_length/2.0)) != (root_shift[2] - (root_scale/2.0))
#endif
		   ) { /*This cell not adjacent to XY plane - make it solid*/
		  child_cells[i] -> cell_type = SOLID;
		  child_cells[i] -> cell_volume = 0;

		  if((directions != ALL) && (child_cells[i] -> cell_type == SOLID)) /*Set flux areas of solid children first*/
		    {
		      for(j = 0; j < NUM_FACES; j++)
			for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
			  child_cells[i] -> flux_area[j][k] = 0;
		    }
		}
	      }
	    }
	    else {
	      if(Parent -> Solid != NULL)
		{
		  for(i = 0; i < MAX_NUM_CHILDREN; i++)
		    child_cells[i] -> Solid = Parent -> Solid;

		  /*If parent only intersects (even potentially so) 1 body, children can only potentially intersect this body.  1/more
		    bodies we'll use the lists to save storage*/
		}          
      	              
	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		{
		  if(child_cells[i] -> cell_type == UNKNOWN)
		    {
		      if(Parent_type == SOLID) /*Parent SOLID => children SOLID*/
			{
			  child_cells[i] -> cell_type = SOLID;
			  child_cells[i] -> cell_volume = 0;
			}
		      else if(Parent_type == FLUID) /*Parent FLUID => children FLUID*/
			{
			  child_cells[i] -> cell_type = FLUID;
			  child_cells[i] -> cell_volume = child_fluid_volume;
			}
		      else /*So parent's type is UNKNOWN/INTERSECTED*/
			{		  
			  if((refine_steps == 1) || (use_geom_engine_always == TRUE))
			    { 
			      if(set_cell_type(child_cells[i], child_fluid_volume, child_edge_length, 
					       &(child_cells[i] -> Intersected_bbox), Parent -> Solid) == ERROR)
				return(ERROR);
			  
			      /*set_cell_type() will also determine if child_cells[i] -> Solid is NULL or not*/
		      
			      if(Parent -> Solid == NULL) /*Then the Intersected_bbox[i] would be built, else the list wouldn't be built
							    as Parent -> Solid would just be used.*/
				child_cells[i] -> bbox_unknown = FALSE;
			    }
			}
		    }

		  if((directions != ALL) && (child_cells[i] -> cell_type == SOLID)) /*Set flux areas of solid children first*/
		    {
		      for(j = 0; j < NUM_FACES; j++)
			for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
			  child_cells[i] -> flux_area[j][k] = 0;
		    }

		  if((child_cells[i] -> cell_type == INTERSECTED) && 
		     (child_cells[i] -> cell_volume < SMALL_CELL_LIMIT*child_fluid_volume)) /*Child cell is small cell*/
		    child_cells[i] -> is_small_cell = TRUE;
		  else child_cells[i] -> is_small_cell = FALSE;
		}
	    }
	  }
      
      
      /*\b Proceed to set the flux areas of children*/

      if(do_area == TRUE)
	{
	  child_face_area = CALC_CELL_FLUID_AREA(child_cells[0]);      

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

	  for(i = 0; i < MAX_NUM_CHILDREN; i++) /*Cycle through each child*/ 
	    {
	      if(child_cells[i] -> cell_type == SOLID) /*If cell type is SOLID, it can't flux at all*/
		{ 
		  for(dir = 0; dir < num_faces; dir++) /*Cycle through each direction*/
		    {
		      if(num_faces == 3)
			j = dirs[dir];
		      else j = dir;

		      /*We should also set the flux area of the corresponding neighbour(s) if they haven't already
			been set*/

                      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) /*Cycle through each quadrant*/
                        child_cells[i] -> flux_area[j][k] = 0;
		      
		      if(child_cells[i] -> face_neighbours[j][0] != NULL)
			push_neighbour_flux_areas(child_cells[i], j, child_edge_length);		
		    }		      
		}
	      else /*The cell is FLUID/INTERSECTED/UNKNOWN*/		    
		{ 
		  /*Before interrogating geometry engine, see if useful info can be gleaned from neighbours*/

		  switch(i) /*Non-sibling directions vary with cell number*/
		    {
		    case 0: /*Cell is child 0*/
		      nsdirection[0] = STH; nsdirection[1] = WST; nsdirection[2] = LWR; break;
		    case 1: 
		      nsdirection[0] = NTH; nsdirection[1] = WST; nsdirection[2] = LWR; break;
		    case 2:
		      nsdirection[0] = STH; nsdirection[1] = EST; nsdirection[2] = LWR; break;
		    case 3:
		      nsdirection[0] = NTH; nsdirection[1] = EST; nsdirection[2] = LWR; break;
		    case 4:
		      nsdirection[0] = STH; nsdirection[1] = WST; nsdirection[2] = UPR; break;
		    case 5:
		      nsdirection[0] = NTH; nsdirection[1] = WST; nsdirection[2] = UPR; break;
		    case 6:
		      nsdirection[0] = STH; nsdirection[1] = EST; nsdirection[2] = UPR; break;
		    case 7:
		      nsdirection[0] = NTH; nsdirection[1] = EST; nsdirection[2] = UPR; break;
		    }
		  
		  for(dir = 0; dir < num_faces; dir++) /*Cycle through each direction*/
		    {
		      if(num_faces == 3)
			j = dirs[dir];
		      else j = dir;

		      switch(j)
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
		      		      
		      if((j == nsdirection[0]) || (j == nsdirection[1]) || (j == nsdirection[2]))
			{ 			      
			  face_set = FALSE;

			  /*There won't exist faces with some quadrants with known flux areas and some unknown as we
			    deal with whole faces at a time.  We 1st check neighbours as we want shared areas to match.*/

			  if(
#if FLOAT_EQ_STABLE
			     EQ(child_cells[i] -> flux_area[j][0], UNKNOWN)
#else
			     child_cells[i] -> flux_area[j][0] == UNKNOWN
#endif
			     )
			    {
			      if(child_cells[i] -> face_neighbours[j][0] != NULL)
				{
				  if(child_cells[i] -> face_neighbours[j][0] -> cell_level == (child_level-1)) 
				    { /*1 level coarser*/

				      switch(j)
					{
					case EST: /*Cells on east face are 2,3,6,7*/ 
					  if(i <= 3)
					    quadrant_num = i-2;
					  else quadrant_num = i-4;
					  break;
					case WST:
					  if(i <= 1)
					    quadrant_num = i;
					  else quadrant_num = i-2;
					  break;
					case NTH:
					  quadrant_num = (i-1)/2; break;
					case STH:
					  quadrant_num = i/2; break;
					case UPR:
					  quadrant_num = i-4; break;
					case LWR:
					  quadrant_num = i; break;
					}

				      if(
#if FLOAT_EQ_STABLE
					 EQ(child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][quadrant_num], 0)
#else
					 child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][quadrant_num] == 0
#endif
					 )
					{
					  for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) 
					    child_cells[i] -> flux_area[j][k] = 0; /*The whole face is obstructed*/

					  face_set = TRUE;

#if AXI_OK
					  if((twoD == 2) && ((j == EST) || (j == WST))) {
					    child_cells[i] -> face_r_centroid[j][0] = 0;
					    child_cells[i] -> face_r_centroid[j][1] = 0;
					  }
#endif
					}
				      else {

					if(twoD != 0) { /*Must use edge lengths instead of 'areas'*/
					  if(
#if FLOAT_EQ_STABLE
					     EQ(child_cells[i]->face_neighbours[j][0]->flux_area[opp_dir][quadrant_num], 
						child_edge_length)
#else
					     child_cells[i]->face_neighbours[j][0]->flux_area[opp_dir][quadrant_num] == 
					     child_edge_length 
#endif
					    ) {
					    child_cells[i] -> flux_area[j][0] = child_edge_length/2.0;
					    child_cells[i] -> flux_area[j][1] = child_edge_length/2.0;
					    child_cells[i] -> flux_area[j][2] = 0;
					    child_cells[i] -> flux_area[j][3] = 0;
				  
					    face_set = TRUE;

#if AXI_OK
					    if((twoD == 2) && ((j == EST) || (j == WST))) {
					      child_cells[i]->face_r_centroid[j][0] = child_cells[i]->centroid[1] - child_edge_length/4.0;
					      child_cells[i]->face_r_centroid[j][1] = child_cells[i]->centroid[1] + child_edge_length/4.0;
					    }
#endif 
					  }
					}
					else {
					  if(
#if FLOAT_EQ_STABLE
					     EQ(child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][quadrant_num],
						child_face_area)
#else
					     child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][quadrant_num] == 
					     child_face_area
#endif
					     )
					    { /*An unobstructed quadrant should have the same area as the corresponding child's face*/

					      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
						child_cells[i] -> flux_area[j][k] = child_face_area/MAX_NUM_NEIGHBOURS;
					      /*All quadrants of the children are thus unobstructed too*/
				  
					      face_set = TRUE;
					    }
					}
				      }
				    }
				  else if(((child_cells[i] -> face_neighbours[j][0] -> cell_level) == child_level) &&
#if FLOAT_EQ_STABLE
					  NOT_EQ(child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][0], UNKNOWN)
#else
					  (child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][0] != UNKNOWN)
#endif
					  )
				    { 
				      /*If actual neighbour is at same level, we might be able to use their flux areas*/
			      
				      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) 
					{
					  child_cells[i] -> flux_area[j][k] = child_cells[i] -> face_neighbours[j][0] ->
					    flux_area[opp_dir][k];				      
					}
				      face_set = TRUE;

#if AXI_OK
				      if((twoD == 2) && ((j == EST) || (j == WST))) {
					child_cells[i] -> face_r_centroid[j][0] = child_cells[i] -> face_neighbours[j][0] ->
					  face_r_centroid[opp_dir][0];
					child_cells[i] -> face_r_centroid[j][1] = child_cells[i] -> face_neighbours[j][0] ->
					  face_r_centroid[opp_dir][1];
				      }
#endif				      
				    }
				  else if(child_cells[i] -> face_neighbours[j][0] -> cell_level == (child_level + 1))
				    { /*If the finer neighbour's quadrant areas are known, we can sum them to form the
					face area of an adjacent child*/

				      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) /*Go through each quadrant of the child*/
					{
					  if(
#if FLOAT_EQ_STABLE
					     NOT_EQ(child_cells[i] -> face_neighbours[j][k] -> flux_area[opp_dir][0], UNKNOWN)
#else
					     child_cells[i] -> face_neighbours[j][k] -> flux_area[opp_dir][0] != UNKNOWN
#endif
					     )
					    {
					      child_cells[i] -> flux_area[j][k] = 
						(child_cells[i] -> face_neighbours[j][k] -> flux_area[opp_dir][0] + 
						 child_cells[i] -> face_neighbours[j][k] -> flux_area[opp_dir][1] + 
						 child_cells[i] -> face_neighbours[j][k] -> flux_area[opp_dir][2] + 
						 child_cells[i] -> face_neighbours[j][k] -> flux_area[opp_dir][3]);
					    }
					}

				      face_set = TRUE;

#if AXI_OK
				      if((twoD == 2) && ((j == EST) || (j == WST))) {
					for(k = 0; k < 2; k++) {
					  if((child_cells[i]->face_neighbours[j][k]->flux_area[opp_dir][0] > 0) ||
					     (child_cells[i]->face_neighbours[j][k]->flux_area[opp_dir][1] > 0))
					    child_cells[i] -> face_r_centroid[j][k] =
					      ((child_cells[i]->face_neighbours[j][k]->flux_area[opp_dir][0])*
					       (child_cells[i]->face_neighbours[j][k]->face_r_centroid[opp_dir][0]) + 
					       (child_cells[i]->face_neighbours[j][k]->flux_area[opp_dir][1])*
					       (child_cells[i]->face_neighbours[j][k]->face_r_centroid[opp_dir][1]))/
					      ((child_cells[i]->face_neighbours[j][k]->flux_area[opp_dir][0])+
					       (child_cells[i]->face_neighbours[j][k]->flux_area[opp_dir][1]));
					  else child_cells[i] -> face_r_centroid[j][k] = 0;
					} 
				      }
#endif
				    }
				}
			    }
			  else face_set = TRUE; /*Flux areas have already been set by another neighbour cell*/

			  if(face_set == FALSE) /*So unable to pull flux area from non-sibling neighbour*/
			    {	
			      /*Might still be able to use parent's flux areas*/	
			      
			      switch(j)
				{
				case NTH: /*Cells on north face are 1,3,5,7*/
				  quadrant_num = (i-1)/2; break; 
				case STH:
				  quadrant_num = i/2; break;
				case EST: /*Cells on east face are 2,3,6,7*/ 
				  if(i <= 3)
				    quadrant_num = i-2;
				  else quadrant_num = i-4;
				  break; 
				case WST:
				  if(i <= 1)
				    quadrant_num = i;
				  else quadrant_num = i-2;
				  break;
				case LWR: 
				  quadrant_num = i; break;
				case UPR: 
				  quadrant_num = i-4; break;
				}
			      
			      if(twoD != 0) {
				/*We always set all XY (or LWR) areas to be 0 of parents for 2D modes, so things must be dealt with
				  differently here*/

				if((Parent->cell_type == FLUID) && (j == LWR)) {  /*I.e. parent has no obstructed XY area*/
				  for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
				    child_cells[i] -> flux_area[j][k] = child_face_area/MAX_NUM_NEIGHBOURS;

				  face_set = TRUE;

#if AXI_OK
				  if(twoD == 2)
				    child_cells[i] -> r_centroid = child_cells[i] -> centroid[1]; 
#endif
				}
				else if((
#if FLOAT_EQ_STABLE
					EQ(Parent -> flux_area[j][quadrant_num], 0)
#else
					Parent -> flux_area[j][quadrant_num] == 0
#endif
					) && (j != LWR)) { /*As LWR flux areas are all 0, but these areas still need to be interrogated*/
				  for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) 
				    child_cells[i] -> flux_area[j][k] = 0;

				  face_set = TRUE;

#if AXI_OK
				  if((twoD == 2) && ((j == EST) || (j == WST))) {
				    child_cells[i] -> face_r_centroid[j][0] = 0;
				    child_cells[i] -> face_r_centroid[j][1] = 0;
				  }
#endif
				}
				else if(
#if FLOAT_EQ_STABLE
					EQ(Parent -> flux_area[j][quadrant_num], child_edge_length)
#else
					Parent -> flux_area[j][quadrant_num] == child_edge_length
#endif
					) {
				  child_cells[i] -> flux_area[j][0] = child_edge_length/2.0;
				  child_cells[i] -> flux_area[j][1] = child_edge_length/2.0;
				  child_cells[i] -> flux_area[j][2] = 0;
				  child_cells[i] -> flux_area[j][3] = 0;
				  
				  /*Dimensionality reduced in 2D - we must have two 'redundant' quadrants to correspond
				    to the quadtree*/
		  
				  face_set = TRUE;

#if AXI_OK
				  if((twoD == 2) && ((j == EST) || (j == WST))) {
				    child_cells[i]->face_r_centroid[j][0] = child_cells[i]->centroid[1] - child_edge_length/4.0;
				    child_cells[i]->face_r_centroid[j][1] = child_cells[i]->centroid[1] + child_edge_length/4.0;
				  }
#endif
				}
			      }
			      else { /*Working in standard 3D mode*/
				if(
#if FLOAT_EQ_STABLE
				   EQ(Parent -> flux_area[j][quadrant_num], 0)
#else
				   Parent -> flux_area[j][quadrant_num] == 0
#endif
				   ) /*So parent's quadrant is fully obstructed*/
				  {
				    for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) 
				      child_cells[i] -> flux_area[j][k] = 0; /*The whole face is obstructed*/
				  
				    face_set = TRUE;
				  }
				else if(
#if FLOAT_EQ_STABLE
					EQ(Parent -> flux_area[j][quadrant_num], child_face_area)
#else
					Parent -> flux_area[j][quadrant_num] == child_face_area
#endif
					)
				  { /*An unobstructed quadrant should have the same area as the corresponding child's face*/
				  
				    for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
				      child_cells[i] -> flux_area[j][k] = child_face_area/MAX_NUM_NEIGHBOURS;
				    /*All quadrants of the children are thus unobstructed too*/
				  
				    face_set = TRUE;
				  }
			      }

			      if((face_set == FALSE) && ((refine_steps == 1) || (use_geom_engine_always == TRUE)))
				{ /*Can't use neighbour and parental values - must now interrogate geometry engine*/

				  if(set_face_flux_area(child_cells[i], j, child_face_area, child_edge_length, 
							&(child_cells[i] -> Intersected_bbox), NULL, 
							child_cells[i] -> bbox_unknown) == ERROR)
				    return(ERROR);
			      			      
				  child_cells[i] -> bbox_unknown = FALSE; 
				  face_set = TRUE;

				  if((twoD != 0) && (j != LWR))  /*Must modify quadrant 'areas' to become line lengths*/
				    for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
				      child_cells[i] -> flux_area[j][k] = (child_cells[i]->flux_area[j][k])/twoD_thickness;
				  /*Quadrants 2 and 3 should be zero if geometry set up correctly*/
				}

			      /*Now set flux areas of neighbours - which will be unknown (else we would've used them)*/

			      if(child_cells[i] -> face_neighbours[j][0] != NULL) 
				push_neighbour_flux_areas(child_cells[i], j, child_edge_length);
			    }
			}
		      else
			{ /*Looking in sibling direction*/  
			  
			  face_set = FALSE;

			  if(
#if FLOAT_EQ_STABLE
			     NOT_EQ(child_cells[i] -> flux_area[j][0], UNKNOWN)
#else
			     child_cells[i] -> flux_area[j][0] != UNKNOWN
#endif
			     )
			    face_set = TRUE;
			    			  
			  if((face_set == FALSE) && 
#if FLOAT_EQ_STABLE
			     NOT_EQ(child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][0], UNKNOWN)
#else
			     (child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][0] != UNKNOWN)
#endif
			     )
			    {				  
			      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
				{
				  child_cells[i] -> flux_area[j][k] = child_cells[i] -> face_neighbours[j][0] ->
				    flux_area[opp_dir][k];				      
				}
			      face_set = TRUE;

#if AXI_OK
			      if((twoD == 2) && ((j == EST) || (j == WST))) {
				child_cells[i] -> face_r_centroid[j][0] = child_cells[i] -> face_neighbours[j][0] ->
				  face_r_centroid[opp_dir][0];
				child_cells[i] -> face_r_centroid[j][1] = child_cells[i] -> face_neighbours[j][0] ->
				  face_r_centroid[opp_dir][1];
			      }
#endif			      
			    }
			  
			  if(face_set == FALSE)
			    {
			      if(Parent_type == FLUID)
				{
				  /*If the parent is FLUID then all children cells can exchange fluxes in sibling directions
				    (assuming no infinitely thin plates).*/
				  
				  if(twoD != 0) {
				    child_cells[i] -> flux_area[j][0] = child_edge_length/2.0;
				    child_cells[i] -> flux_area[j][1] = child_edge_length/2.0;
				    child_cells[i] -> flux_area[j][2] = 0;
				    child_cells[i] -> flux_area[j][3] = 0;
				  }
				  else {
				    for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) 
				      child_cells[i] -> flux_area[j][k] = child_face_area/MAX_NUM_NEIGHBOURS;
				  }
			      
				  face_set = TRUE;

#if AXI_OK
				  if((twoD == 2) && ((j == EST) || (j == WST))) {
				    child_cells[i]->face_r_centroid[j][0] = child_cells[i]->centroid[1] - child_edge_length/4.0;
				    child_cells[i]->face_r_centroid[j][1] = child_cells[i]->centroid[1] + child_edge_length/4.0;
				  }
#endif
				}

			      if((face_set == FALSE) && ((refine_steps == 1) || (use_geom_engine_always == TRUE)))
				{
				  /*Must interrogate geometry engine again*/
			      
				  if(Parent -> Solid != NULL)
				    {
				      if(set_face_flux_area(child_cells[i], j, child_face_area, child_edge_length, 
							    &(child_cells[i] -> Intersected_bbox), Parent -> Solid,  FALSE) == ERROR)
					return(ERROR);
				    }
				  else 
				    {
				      if(set_face_flux_area(child_cells[i], j, child_face_area, child_edge_length,
							    &(child_cells[i] -> Intersected_bbox), NULL, 
							    child_cells[i] -> bbox_unknown) == ERROR)
					return(ERROR); 
				      child_cells[i] -> bbox_unknown = FALSE;

				      if(twoD != 0)  /*Must modify quadrant 'areas' to become line lengths*/
					for(k = 0; k < MAX_NUM_NEIGHBOURS; k++)
					  child_cells[i] -> flux_area[j][k] = (child_cells[i]->flux_area[j][k])/twoD_thickness;
				      /*Quadrants 2 and 3 should be zero if geometry set up correctly*/
				    }
			      
				  /*If the parent only intersects 1 body, then sibling fluxes can only be obstructed by this body*/
			      
				  face_set = TRUE;
				}

			      /*We should update corresp. neighbour if it has UNKNOWN flux areas too (which it will
				or else we would've used it earlier)*/

			      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) 
				{
				  child_cells[i] -> face_neighbours[j][0] -> flux_area[opp_dir][k] = 
				    child_cells[i] -> flux_area[j][k];  
				}

#if AXI_OK
			      if((twoD == 2) && ((j == EST) || (j == WST))) {
				child_cells[i] -> face_neighbours[j][0] -> face_r_centroid[opp_dir][0] = 
				  child_cells[i] -> face_r_centroid[j][0];
				child_cells[i] -> face_neighbours[j][0] -> face_r_centroid[opp_dir][1] = 
				  child_cells[i] -> face_r_centroid[j][1];
			      }
#endif
			    }
			}
		    }

#if 0
		  if(wall_rep == 1) /*For smoothed wall representations, we need to compute total obstructed area and surface normal*/
		    {
		      for(j = 0; j < 3; j++)
			{
			  k = j*2; /*'Positive' directions along the 3 axes*/
			  opp_dir = k+1; /*'Negative' directions along the 3 axes*/

			  child_cells[i]->wall_norm[j] = 
			    (child_cells[i]->flux_area[opp_dir][0]+child_cells[i]->flux_area[opp_dir][1]+
			     child_cells[i]->flux_area[opp_dir][2]+child_cells[i]->flux_area[opp_dir][3]) - 
			    (child_cells[i]->flux_area[k][0]+child_cells[i]->flux_area[k][1]+
			     child_cells[i]->flux_area[k][2]+child_cells[i]->flux_area[k][3]);
			}

		      child_cells[i] -> wall_area = SQRT(SQR(child_cells[i]->wall_norm[0]) + SQR(child_cells[i]->wall_norm[1]) +
							 SQR(child_cells[i]->wall_norm[2])); /*Area is magnitude of unnormalized normal*/
		      for(j = 0; j < 3; j++)
			child_cells[i] -> wall_norm[j] = (child_cells[i]->wall_norm[j])/(child_cells[i]->wall_area); /*Normalize*/
		    }
#endif

		} /*Face areas worth setting*/
	    } /*End cycle through children*/	  

	  if((twoD != 0) && ((directions == ALL) || (directions == NEGATIVE))) { 
	    /*Must reduce dimensionality by making 'volume' = xy area, and 'area' = line length, i.e. volume per
	      unit thickness and area per unit thickness*/

	    for(i = 0; i < MAX_NUM_CHILDREN; i++) {
	      if(child_cells[i] -> cell_type != SOLID) { /*'Upper' children will be always 'solid' type*/

		child_cells[i] -> cell_volume = (child_cells[i]->flux_area[LWR][0] + child_cells[i]->flux_area[LWR][1] + 
						 child_cells[i]->flux_area[LWR][2] + child_cells[i]->flux_area[LWR][3]);

		if(
#if FLOAT_EQ_STABLE
		   EQ(child_cells[i]->cell_volume, child_face_area)
#else
		   child_cells[i] -> cell_volume == child_face_area
#endif
		   ) {
		  child_cells[i] -> cell_type = FLUID;
		}
		else if(
#if FLOAT_EQ_STABLE
			EQ(child_cells[i]->cell_volume, 0)
#else
			child_cells[i] -> cell_volume == 0
#endif
			) {
		  child_cells[i] -> cell_type = SOLID;

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

		  for(dir = 0; dir < num_faces; dir++) /*Cycle through each direction*/
		    {
		      if(num_faces == 3)
			j = dirs[dir];
		      else j = dir;

		      /*We should reset corresponding neighbours' areas*/

		      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) /*Cycle through each quadrant*/
			child_cells[i] -> flux_area[j][k] = 0;

		      if(child_cells[i] -> face_neighbours[j][0] != NULL)
			push_neighbour_flux_areas(child_cells[i], j, child_edge_length);		
		    }	      
		}
		else {
		  child_cells[i] -> cell_type = INTERSECTED;	      
		}

		if(child_cells[i] -> cell_type != SOLID) {

		  if(child_cells[i] -> cell_volume < (SMALL_CELL_LIMIT*child_face_area)) /*Child cell is small cell*/
		    child_cells[i] -> is_small_cell = TRUE;
		  else child_cells[i] -> is_small_cell = FALSE;

		  /*UPR face areas should already be 0.  Now reset LWR face areas to be 0*/
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		    child_cells[i]->flux_area[LWR][j] = 0;
		}
	      }
	    }
	  }

	  /*End of flux area setting section*/
	}            
      

      /*\b Proceed to set additional_refine_times of the children (TRUE/FALSE)*/

      if(do_coarsenability == TRUE)
	{
	  short int Parent_additional_refine_times = Parent -> additional_refine_times; 
      
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if(child_cells[i] -> additional_refine_times == UNKNOWN) /*Not already set before*/
		{
		  if((child_cells[i] -> cell_type == SOLID) || (child_cells[i] -> cell_type == FLUID))
		    {
		      /*Solid/fluid cells need not be additionally refined for any geometrical reasons unless they're
			too coarse for the initial grid*/
		  
		      if(child_level < min_refinement_level)
			{ 
			  if(allow_multiple_refinements == FALSE)
			    child_cells[i] -> additional_refine_times = 1; /*Can only refine once each time*/
			  else child_cells[i] -> additional_refine_times = min_refinement_level - child_level; 
			  /*Refine until min_refinement_level attained*/
			}
		      else child_cells[i] -> additional_refine_times = 0;
		    }
		  else if(Parent_additional_refine_times == 0)
		    child_cells[i] -> additional_refine_times = 0;
#if NOTINCLUDE
		  else if(Parent_additional_refine_times != UNKNOWN) /*Check parent's data*/ 
		    {
		      if(Parent_additional_refine_times == 0)
			child_cells[i] -> additional_refine_times = 0;
		      else child_cells[i] -> additional_refine_times = (Parent -> additional_refine_times) - 1;
		    }
#endif
		  else 
		    {		      
		      if((use_geom_engine_always == TRUE) || (refine_steps == 1))
			{
			  if(set_cell_additional_refine_times(child_cells[i], &(child_cells[i] -> Intersected_bbox), 
							      child_cells[i] -> Solid) == ERROR)
			    return(ERROR);
			}
		    }
		} /*Child's additional_refine_times is already known*/
	    }
	}


      if(rm_temp_lists == TRUE)
	{
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    destroy_list_bbox(&(child_cells[i] -> Intersected_bbox)); /*No need for these lists anymore*/
      
	  /**< \b NOTE - If I don't deallocate lists, overall CPU performance better - but do I have enough memory?*/
	}
      
      return(NOERROR);
      
    } /*Setting cell type, additional_refine_times and flux areas bypassed if child_level < level_before_can_interrogate*/
}

/*------------------------------------------------------------------*/

/**\brief Set parent's geometrical data (based on children)*/

short int set_parent_geom_data(Cart_cell child_cells[])
{
  short int i, j, k, x;
  Cart_cell Parent = child_cells[0] -> parent;
  double child_fluid_volume = CALC_CELL_FLUID_VOL(child_cells[0]);
  short int Parent_level = Parent -> cell_level;
  short int child_level = child_cells[0] -> cell_level;

  /*Proceed to set cell type of parent (use more accurate children's values)*/
  k = 0;
  Parent -> cell_volume = 0;
  for(i = 0; i < MAX_NUM_CHILDREN; i++)
    Parent -> cell_volume += (child_cells[i] -> cell_volume); /*For 2D cases 'volume' is really 'area'*/

  if(twoD == 0) { /*Min refine level == min interrogate level, so need not worry that a parent will have unknown is_small_cell*/
    if(Parent -> cell_type == UNKNOWN) 
      {
	if(Parent -> cell_volume == 0) /*Hope that exact floating point comparison OK here - bcos subcell volumes very tiny*/
	  Parent -> cell_type = SOLID;
	else if(
#if FLOAT_EQ_STABLE
		EQ(Parent -> cell_volume, MAX_NUM_CHILDREN*child_fluid_volume)
#else
		Parent -> cell_volume == (MAX_NUM_CHILDREN*child_fluid_volume)
#endif
		)
	  Parent -> cell_type = FLUID;
	else Parent -> cell_type = INTERSECTED;

	/*Compare children -> Solid pointers to see if the parent intersects 1 body only*/
      
	if(Parent -> cell_type == FLUID)
	  Parent -> Solid = NULL;
	else 
	  {
	    if(Parent -> Solid == NULL) /*Parent -> Solid not previously set already*/
	      {
		Body tempbody = NULL; 
		j = TRUE;
	      
		for(i = 0; i < MAX_NUM_CHILDREN; i++)
		  {
		    if(child_cells[i] -> cell_type != FLUID) /*Disregard FLUID children*/
		      {
			if(child_cells[i] -> Solid == NULL) /*An intersected/solid child has NULL Solid pointer*/
			  {
			    Parent -> Solid = NULL; /*More than 1 body intersects parent too*/
			    j = FALSE;
			    break;
			  }
			else if(child_cells[i] -> Solid != NULL) /*Only 1 body is stored*/
			  {
			    if(tempbody == NULL) /*Let tempbody be the first Solid body encountered*/
			      tempbody = child_cells[i] -> Solid;
			    else if(tempbody != (child_cells[i] -> Solid)) /*Compare bodies*/
			      {
				Parent -> Solid = NULL;
				j = FALSE;
				break;
			      }
			  }
		      }
		  }
		if(j == TRUE) /*All children cells intersect the same Solid*/
		  Parent -> Solid = tempbody;
	      }
	  }
      }

    if(Parent -> is_small_cell == UNKNOWN)
      {
	if((Parent -> cell_type == INTERSECTED) && (Parent -> cell_volume < SMALL_CELL_LIMIT*child_fluid_volume*MAX_NUM_CHILDREN))
	  Parent -> is_small_cell = TRUE;
	else Parent -> is_small_cell = FALSE;
      }
  }  

  /*Set additional_refine_times of parent*/ 
  
  if(Parent -> additional_refine_times == UNKNOWN)
    {
      if((Parent -> cell_type == SOLID) || (Parent -> cell_type ==  FLUID)) 
	{		  
	  if(Parent_level < min_refinement_level)
	    {
	      if(allow_multiple_refinements == FALSE)
		Parent -> additional_refine_times = 1; 
	      else Parent -> additional_refine_times = min_refinement_level - child_level; 
	    }
	  else Parent -> additional_refine_times = 0; 
	}
      else /*The parent is INTERSECTED - see if children know anything*/
	{
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if(child_cells[i] -> additional_refine_times > 0)
		{
		  Parent -> additional_refine_times = 1;
		  break;
		}

	      /*Don't need to know precise value of additional_refine_times; it's enough that parent can't ever be coarsened*/
	    }
	  
	  if(i == MAX_NUM_CHILDREN) 
	    { 
	      List_bbox Bbox_list = NULL;
	      
	      if(set_cell_additional_refine_times(Parent, &Bbox_list, Parent -> Solid) == ERROR)
		return(ERROR);

	      destroy_list_bbox(&Bbox_list); 
	    }
	}

      if(Parent -> additional_refine_times > 0)
	return(NOERROR); /*It's enough that we can't coarsen this cell*/
    }      
 
  /*Flux areas have to be updated each time as previously fluxable areas may no longer be due to solid degeneracies*/
	      
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

#if AXI_OK
  if(twoD == 2) {
    Parent -> r_centroid = ((Parent -> children[0] -> cell_volume)*(Parent -> children[0] -> r_centroid) + 
			    (Parent -> children[1] -> cell_volume)*(Parent -> children[1] -> r_centroid) + 
			    (Parent -> children[2] -> cell_volume)*(Parent -> children[2] -> r_centroid) +
			    (Parent -> children[3] -> cell_volume)*(Parent -> children[3] -> r_centroid))/
      ((Parent -> children[0] -> cell_volume) + (Parent -> children[1] -> cell_volume) + (Parent -> children[2] -> cell_volume) +
       (Parent -> children[3] -> cell_volume));

    if((Parent -> children[2] -> flux_area[EST][0] > 0) || (Parent -> children[2] -> flux_area[EST][1] > 0))
      Parent -> face_r_centroid[EST][0] = ((Parent->children[2] -> face_r_centroid[EST][0])*(Parent -> children[2] -> flux_area[EST][0]) +
					   (Parent->children[2] -> face_r_centroid[EST][1])*(Parent -> children[2] -> flux_area[EST][1]))/
	((Parent -> children[2] -> flux_area[EST][0]) + (Parent -> children[2] -> flux_area[EST][1]));
    else Parent -> face_r_centroid[EST][0] = 0;

    if((Parent -> children[3] -> flux_area[EST][0] > 0) || (Parent -> children[3] -> flux_area[EST][1] > 0))
      Parent -> face_r_centroid[EST][1] = ((Parent->children[3] -> face_r_centroid[EST][0])*(Parent -> children[3] -> flux_area[EST][0]) +
					   (Parent->children[3] -> face_r_centroid[EST][1])*(Parent -> children[3] -> flux_area[EST][1]))/
	((Parent -> children[3] -> flux_area[EST][0]) + (Parent -> children[3] -> flux_area[EST][1]));
    else Parent -> face_r_centroid[EST][1] = 0;

    if((Parent -> children[0] -> flux_area[WST][0] > 0) || (Parent -> children[0] -> flux_area[WST][1] > 0))
      Parent -> face_r_centroid[WST][0] = ((Parent->children[0] -> face_r_centroid[WST][0])*(Parent -> children[0] -> flux_area[WST][0]) +
					   (Parent->children[0] -> face_r_centroid[WST][1])*(Parent -> children[0] -> flux_area[WST][1]))/
	((Parent -> children[0] -> flux_area[WST][0]) + (Parent -> children[0] -> flux_area[WST][1]));
    else Parent -> face_r_centroid[WST][0] = 0;

    if((Parent -> children[1] -> flux_area[WST][0] > 0) || (Parent -> children[1] -> flux_area[WST][1] > 0))
      Parent -> face_r_centroid[WST][1] = ((Parent->children[1] -> face_r_centroid[WST][0])*(Parent -> children[1] -> flux_area[WST][0]) +
					   (Parent->children[1] -> face_r_centroid[WST][1])*(Parent -> children[1] -> flux_area[WST][1]))/
	((Parent -> children[1] -> flux_area[WST][0]) + (Parent -> children[1] -> flux_area[WST][1]));
    else Parent -> face_r_centroid[WST][1] = 0;
  }
#endif


#if 0
  if(wall_rep == 1)
    {
      for(j = 0; j < 3; j++)
	{
	  k = j*2;
	  x = k+1;

	  Parent -> wall_norm[j] = 
	    (Parent->flux_area[x][0]+Parent->flux_area[x][1]+Parent->flux_area[x][2]+Parent->flux_area[x][3]) -
	    (Parent->flux_area[k][0]+Parent->flux_area[k][1]+Parent->flux_area[k][2]+Parent->flux_area[k][3]);
	}

      Parent -> wall_area = SQRT(SQR(Parent->wall_norm[0]) + SQR(Parent->wall_norm[1]) + SQR(Parent->wall_norm[2]));
      for(j = 0; j < 3; j++)
	Parent->wall_norm[j] = (Parent->wall_norm[j])/(Parent->wall_area);
    }
#endif

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Set cell type and cell volume of cell*/

short int set_cell_type(Cart_cell C, double cell_volume, double edge_length, List_bbox *Intersected_bbox, Body Solid) 
{  
  int i;
  double cell_bbox[2][3], rootub[3], rootlb[3];
  
  /*Get cell bounding box from verticies 0 and 7*/
  
  for(i = 0; i < 3; i++)
    {
#if !VTXM
      cell_bbox[0][i] = C -> verticies[0] -> loc[i];
      cell_bbox[1][i] = C -> verticies[7] -> loc[i];
#else
      cell_bbox[0][i] = C -> centroid[i] - 0.5*CALC_CELL_EDGE_LENGTH(C);
      cell_bbox[1][i] = C -> centroid[i] + 0.5*CALC_CELL_EDGE_LENGTH(C);
#endif      

      rootlb[i] = root_shift[i] - 0.5*root_scale;
      rootub[i] = root_shift[i] + 0.5*root_scale;
    }
    
  if(((*Intersected_bbox) == NULL) && (Solid == NULL))
    {           
      double lower_bound[6]; double upper_bound[6]; /*Bounding box of cell in hyperspace*/

      for(i = 0; i < 6; i++)
	{
	  if(i >= 3)
	    {
	      lower_bound[i] = cell_bbox[0][i%3];
	      upper_bound[i] = rootub[i%3];
	    }
	  else 
	    {
	      lower_bound[i] = rootlb[i];
	      upper_bound[i] = cell_bbox[1][i];
	    }	  
	}
      
      if(traverse_body_ADT(lower_bound, upper_bound, ADTbody, &(*Intersected_bbox)) == ERROR)
	return(ERROR);
    }

  if(((*Intersected_bbox) == NULL) && (Solid == NULL)) /*So cell intersects no bounding boxes - will be completely fluid*/
    {
      C -> cell_type = FLUID;
      C -> cell_volume = cell_volume;
      C -> Solid = NULL;
      return(NOERROR);
    }
  else
    {
      Body temp_body = NULL;
      short int is_unique_body = TRUE;
      double vol, subcell_vol;      

      int num_subcells_inside = 0; /*No. subcells inside intersected bodies*/
      
      traverse_volume_ADT(cell_bbox[0], edge_length, &num_subcells_inside, ADT3Dpoint, 
			  &(*Intersected_bbox), Solid, &temp_body, &is_unique_body, NULL);
                  
      if(num_subcells_inside == 0) 
	{
	  C -> cell_type = FLUID;
	  C -> cell_volume = cell_volume;
	  C -> Solid = NULL;
	}
      else if(num_subcells_inside == CUBE(volume_subcell_num)) /*All subcells immersed*/
	{  
	  C -> cell_type = SOLID;
	  C -> cell_volume = 0.0;
	  
	  if(is_unique_body == TRUE) /*So cell only inside 1 body*/
	    {
	      if(Solid != NULL)
		C -> Solid = Solid; 
	      else C -> Solid = temp_body;
	    }
	  else C -> Solid = NULL; /*Cell is inside more than 1 body*/
	}
      else /**\b NOTE - for the moment if even 1 subcell is inside a body, it's an intersected cell*/
	{
	  C -> cell_type = INTERSECTED;
	  
	  vol = 0.0;
	  subcell_vol = cell_volume/((double) (CUBE(volume_subcell_num)));
	  
	  for(i = 0; i < num_subcells_inside; i++)
	    vol = vol + subcell_vol;
	  
	  C -> cell_volume = cell_volume - vol;
	  
	  if(is_unique_body == TRUE)
	    {
	      if(Solid != NULL)
		C -> Solid = Solid;
	      else C -> Solid = temp_body;
	    }
	  else C -> Solid = NULL; 
	}
      
      return(NOERROR);
    }      
}

/*------------------------------------------------------------------*/

/**\brief Set additional_refine_times of a cell - typically needs curvatures but as for the moment only volume
 representations are used, is a very simple function*/

short int set_cell_additional_refine_times(Cart_cell C, List_bbox *Intersected_bbox, Body Solid) 
{ 
  short int i;

  /*Only INTERSECTED cells are submitted to this routine; SOLID/FLUID cells won't even execute this function - later on
    may need to use Intersected_bbox and Solid; for now not*/

  short int level = C -> cell_level;

  short int n[3] = {0, 0, 0}; short int maxn; 
  
  if(level < min_refinement_level)
    { /*Aat least refine until min_refinement_level reached*/
      
      if(allow_multiple_refinements == FALSE)
   	n[0] = 1;
      else n[0] = min_refinement_level - level;
    }
  
  if(level < min_intersect_refine_level)
    { /*So we should at least refine until min_intersect_refine_level reached*/
      
      if(allow_multiple_refinements == FALSE)
   	n[1] = 1;
      else n[1] = min_intersect_refine_level - level;
    }
  
  if((allow_multiple_refinements == FALSE) && ((n[0] == 1) || (n[1] == 1)))
    {
      C -> additional_refine_times = 1;
      return(NOERROR);
    }     
  else /*So additional_refine_times is UNKNOWN or allow multiple refinements*/
    {
      /* n[2] = max(N), N is no. times cell should be refined given curvatures of Nth body

	If n[2] > 0, include the code fragment below - 

	if(n[2] > 0) {	
	 if(allow_multiple_refinements == FALSE) {
  	  n[2] = 1;
	 }
	}
	
	For the moment we don't have any reason to set n[2] > 0.  n[2] shouldn't also exceed max_refinement_level.*/
    }

  maxn = n[0];
  for(i = 0; i < 3; i++)
    {
      if(n[i] > maxn)
	maxn = n[i];
    }
  
  C -> additional_refine_times = maxn;

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Set the flux areas for all quadrants on a specified face of a cell*/

short int set_face_flux_area(Cart_cell C, short int direction, double unobstructed_face_area, double edge_length, 
			     List_bbox *Intersected_bbox, Body Solid, short int check_intersection) 
{ 
  /*Direction is the direction (N,S,E,W,L,U) corresponding to face*/

  int i;      
  
  if((check_intersection == TRUE) && (Solid == NULL)) 
    {
      double cell_bbox[2][3], rootub[3], rootlb[3];
      double lower_bound[6]; double upper_bound[6]; /*Bounding box of cell in hyperspace*/
      
      for(i = 0; i < 3; i++)
	{
#if !VTXM
	  cell_bbox[0][i] = C -> verticies[0] -> loc[i];
	  cell_bbox[1][i] = C -> verticies[7] -> loc[i];
#else
	  cell_bbox[0][i] = C -> centroid[i] - 0.5*CALC_CELL_EDGE_LENGTH(C);
	  cell_bbox[1][i] = C -> centroid[i] + 0.5*CALC_CELL_EDGE_LENGTH(C);
#endif

	  rootlb[i] = root_shift[i] - 0.5*root_scale;
	  rootub[i] = root_shift[i] + 0.5*root_scale;
	}
                  
      for(i = 0; i < 6; i++)
	{
	  if(i >= 3)
	    {
	      lower_bound[i] = cell_bbox[0][i%3];
	      upper_bound[i] = rootub[i%3];
	    }
	  else 
	    {
	      lower_bound[i] = rootlb[i];
	      upper_bound[i] = cell_bbox[1][i];	  
	    }	  
	}

      if(traverse_body_ADT(lower_bound, upper_bound, ADTbody, &(*Intersected_bbox)) == ERROR)
	return(ERROR);
    }

  if(((*Intersected_bbox) == NULL) && (Solid == NULL)) /*Face intersects no potential bodies*/
    {
      for(i = 0; i < MAX_NUM_NEIGHBOURS; i++) /*All quadrants completely unobstructed*/
	C -> flux_area[direction][i] = unobstructed_face_area/MAX_NUM_NEIGHBOURS;
      return(NOERROR);
    }
  else
    {
      /*Determine if each candidate even intersects face plane in the direction normal to it*/
      
      double face_bbox[3];
      
      short int j;
      int num_subcells_inside[4] = {0, 0, 0, 0}; /*Counter for no. of subcells (in each quadrant) inside intersected bodies*/
      int num2D_subcells_inside[2] = {0, 0};
      double qarea;
      double subcell_area = unobstructed_face_area/((double) (SQR(area_subcell_num)));
      short int axis; /*Identifies which axis is the face plane normal to*/
      double corner_bbox2D[2]; /*Lower corner of face bounding box in 2D*/

      double r_solid[4] = {0, 0, 0, 0};
      double r_avg, r_basis, A_basis, A_free, A_ob;
      
#if !VTXM
      switch(direction) /*Extract relevant co-ordinates of verticies to form lower corner of face's bounding box in 3D*/
	{
	case NTH:
	  v = 1; break;
	case STH:
	  v = 0; break;
	case EST:
	  v = 2; break;
	case WST:
	  v = 0; break;
	case LWR:
	  v = 0; break;
	case UPR:
	  v = 4; break;
	}
      
      for(i = 0; i < 3; i++)
	face_bbox[i] = C -> verticies[v] -> loc[i];
#else
      
      short int nt[3];

      switch(direction) 
	{
	case NTH:
	  nt[0] = 0; nt[1] = 1; nt[2] = 0; break;
	case STH:
	  nt[0] = 0; nt[1] = 0; nt[2] = 0; break;
	case EST:
	  nt[0] = 1; nt[1] = 0; nt[2] = 0; break;
	case WST:
	  nt[0] = 0; nt[1] = 0; nt[2] = 0; break;
	case LWR:
	  nt[0] = 0; nt[1] = 0; nt[2] = 0; break;
	case UPR:
	  nt[0] = 0; nt[1] = 0; nt[2] = 1; break;
	}
   	   
      for(i = 0; i < 3; i++)
	face_bbox[i] = (C -> centroid[i]) + (2*nt[i]-1)*(C -> cell_length)/2.0;
      
#endif
      A_basis = 0;
      if((direction == EST) || (direction == WST))
	{
	  axis = 0;
	  corner_bbox2D[0] = face_bbox[1];
	  corner_bbox2D[1] = face_bbox[2];
	}
      else if((direction == NTH) || (direction == STH))
	{
	  axis = 1;
	  corner_bbox2D[0] = face_bbox[0];
	  corner_bbox2D[1] = face_bbox[2];
	}
      else 
	{
	  axis = 2;
	  corner_bbox2D[0] = face_bbox[0];
	  corner_bbox2D[1] = face_bbox[1];
	}
      
      traverse_area_ADT(axis, face_bbox, corner_bbox2D, face_bbox[axis], edge_length, num_subcells_inside, 
			ADT2Dpoint, &(*Intersected_bbox), Solid, r_solid, num2D_subcells_inside);
      
      for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
	{
	  if(num_subcells_inside[j] == 0) /*Bounding boxes are intersected, but no subcell inside any body*/
	    C -> flux_area[direction][j] = unobstructed_face_area/MAX_NUM_NEIGHBOURS; /*Quadrant completely unobstructed*/
	  else if(num_subcells_inside[j] == ((int) (SQR(area_subcell_num)))/MAX_NUM_NEIGHBOURS) /*All subcells immersed*/
	    C -> flux_area[direction][j] = 0; /*Quadrant completely obstructed*/
	  else /**\b NOTE - for the moment even if 1 subcell is obstructed, quadrant is partially obstructed*/
	    {
	      qarea = 0.0;		  
	      
	      for(i = 0; i < num_subcells_inside[j]; i++)
		qarea = qarea + subcell_area;
	      
	      C -> flux_area[direction][j] = (unobstructed_face_area/MAX_NUM_NEIGHBOURS) - qarea;
	    }
	}

#if AXI_OK
      if(twoD == 2) {
	if((direction == EST) || (direction == WST)) { /*Need to obtain average interface r*/

	  if((num2D_subcells_inside[0] > 0) || (num2D_subcells_inside[1] > 0)) 
	    A_basis = 0.5*(C -> cell_length)*twoD_thickness; /*Quadrant 'basis area' (in 2D)*/
	  
	  if(num2D_subcells_inside[0] > 0) { /*Obstruction in quadrant 0*/
	    if(
#if FLOAT_EQ_STABLE
	       EQ(C -> flux_area[direction][0], 0)
#else
	       C -> flux_area[direction][0] == 0
#endif
	       )
	      C -> face_r_centroid[direction][0] = 0;
	    else {
	      r_basis = (C -> centroid[1]) - 0.25*(C -> cell_length); /*Basis average r*/
	      r_avg = (r_solid[0])/num2D_subcells_inside[0]; /*Obstructed average r*/
	      A_ob = A_basis - (C -> flux_area[direction][0]); /*Obstructed area*/
	      C -> face_r_centroid[direction][0] = (r_basis*A_basis - r_avg*A_ob)/(C -> flux_area[direction][0]);
	    }
	  }
	  else C -> face_r_centroid[direction][0] = (C -> centroid[1]) - 0.25*(C -> cell_length);

	  if(num2D_subcells_inside[1] > 0) {
	    if(
#if FLOAT_EQ_STABLE
	       EQ(C -> flux_area[direction][1], 0)
#else
	       C -> flux_area[direction][1] == 0
#endif
	       )
	      C -> face_r_centroid[direction][1] = 0;
	    else {
	      r_basis = (C -> centroid[1]) + 0.25*(C -> cell_length);
	      r_avg = (r_solid[1])/num2D_subcells_inside[1]; 
	      A_ob = A_basis - (C -> flux_area[direction][1]); 
	      C -> face_r_centroid[direction][1] = (r_basis*A_basis - r_avg*A_ob)/(C -> flux_area[direction][1]);
	    }
	  }
	  else C -> face_r_centroid[direction][1] = (C -> centroid[1]) + 0.25*(C -> cell_length);
	}
	else if(direction == LWR) {

	  if((num_subcells_inside[0] > 0) || (num_subcells_inside[1] > 0) || (num_subcells_inside[2] > 0) ||
	     (num_subcells_inside[3] > 0)) {
	    A_free = (C -> flux_area[LWR][0]) + (C -> flux_area[LWR][1]) + (C -> flux_area[LWR][2]) + (C -> flux_area[LWR][3]);

	    if(
#if FLOAT_EQ_STABLE
	       EQ(A_free, 0)
#else
	       A_free == 0
#endif
	       )
	      C -> r_centroid = 0;
	    else {
	      r_avg = (r_solid[0] + r_solid[1] + r_solid[2] + r_solid[3])/
		(num_subcells_inside[0] + num_subcells_inside[1] + num_subcells_inside[2] + num_subcells_inside[3]);
	      A_ob = unobstructed_face_area - A_free;
	      C -> r_centroid = ((C -> centroid[1])*unobstructed_face_area - r_avg*A_ob)/A_free;
	    }	      
	  }
	  else C -> r_centroid = C -> centroid[1];
	}
      }
#endif
      
      return(NOERROR);
    }
}

/*------------------------------------------------------------------*/

/**\brief Sets IC_flags of cell depending on their relation to any special IC volumes*/

short int set_IC_flag(Cart_cell child_cells[])
{
  int i, j, k, num_subcells_inside, tot_subcells;
  Cart_cell Parent = child_cells[0] -> parent;
  double cell_bbox[2][3];
  double child_edge_length;
  double bbox[6];
  
  List_bbox dummy_list;
  Body dummy_body;
  short int IC_flag = -1;

  if(num_ICs == 0)
    return(NOERROR);
  
  if(child_cells[0] -> cell_level < level_before_can_interrogate) /*Won't bother if cell is too coarse*/ 
    return(NOERROR);  

  tot_subcells = CUBE(volume_subcell_num); /*To test if cell intersected by IC volume, compare how many subcells inside*/
  child_edge_length = CALC_CELL_EDGE_LENGTH(child_cells[0]);

  for(i = 0; i < MAX_NUM_CHILDREN; i++) {
    if(
#if FLOAT_EQ_STABLE
       EQ(child_cells[i] -> IC_flag, UNKNOWN)
#else
       child_cells[i] -> IC_flag == UNKNOWN
#endif
       ) {
      if(Parent -> cell_type == SOLID)
	child_cells[i] -> IC_flag = IRRELEVANT;
      else if(
#if FLOAT_EQ_STABLE
	      EQ(Parent -> IC_flag, OUTSIDE)
#else
	      Parent -> IC_flag == OUTSIDE
#endif
	      )
	child_cells[i] -> IC_flag = OUTSIDE;
      else if((Parent -> IC_flag >= 0) &&
#if FLOAT_EQ_STABLE
	      EQ(ROUND(Parent -> IC_flag), (Parent -> IC_flag))
#else 
	      (ROUND(Parent -> IC_flag) == (Parent -> IC_flag))
#endif
	      ) 
	child_cells[i] -> IC_flag = Parent -> IC_flag; /*If parent's IC_flag is integer, will be wholly inside a volume*/
      else { /*Parent knows nothing, or is intersected by an IC volume (a cell would only see the first IC volume it's
	       intersected by - but children can intersect more)*/
	    
	for(j = 0; j < 3; j++) {
	  cell_bbox[0][j] = child_cells[i] -> verticies[0] -> loc[j];
	  cell_bbox[1][j] = child_cells[i] -> verticies[7] -> loc[j];
	}

	child_cells[i] -> IC_flag = OUTSIDE;

	for(j = 0; j < num_ICs; j++) { /*Check if cell in any one of the special ICs*/
	  for(k = 0; k < 3; k++)
	    if((cell_bbox[0][k] >= IC_bboxes[j][1][k]) || (cell_bbox[1][k] <= IC_bboxes[j][0][k]))
	      break;

	  if(k == 3) { /*Bounding box intersection exists - now test for actual intersection*/
	    num_subcells_inside = 0;

	    for(k = 0; k < 3; k++) {
	      bbox[k] = IC_bboxes[j][0][k];
	      bbox[k+3] = IC_bboxes[j][1][k];
	    }

	    traverse_volume_ADT(cell_bbox[0], child_edge_length, &num_subcells_inside, ADT3Dpoint, &dummy_list,
				IC_regions[j], &dummy_body, &IC_flag, bbox);
	    
	    if(num_subcells_inside > 0) {
	      
#if DETON
	      if(Deton.num_deton_pts > 0)
		child_cells[i] -> IC_flag = (double) j + 0.5; /*If finite detonation is simulated, we should use highest resolution
								anyway*/
	      else {
		if(num_subcells_inside < tot_subcells)
		  child_cells[i] -> IC_flag = (double) j + 0.5; /*Non-integers specified intersected by jth (integer part) IC region*/
		else {
#if FINE_IC
		  child_cells[i] -> IC_flag = (double) j + 0.5;
#else
		  child_cells[i] -> IC_flag = (double) j;
#endif
	}
	      }
#else
	      if(num_subcells_inside < tot_subcells)
		child_cells[i] -> IC_flag = (double) j + 0.5; /*Non-integers specified intersected by jth (integer part) IC region*/
	      else child_cells[i] -> IC_flag = (double) j;
#endif
	      
	      break; /*Inside/intersected by the first seen volume; it is enough, exit*/
	    }
	  }
	}
      }
    }
  }    

  return(NOERROR);  
}

/*------------------------------------------------------------------*/
/**\brief 'Push', or set a neighbour's flux areas in a given direction given the cell's own flux areas are known.
   Assumes a neighbour already exists*/

void push_neighbour_flux_areas(Cart_cell C, short int dir, double child_edge_length)
{
  short int i = C -> child_num;
  short int k, opp_dir;
  double child_face_area = SQR(child_edge_length);
  k = 0; opp_dir = 0;
  switch(dir)
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

  if(C -> face_neighbours[dir][0] -> cell_level == (C -> cell_level - 1)) {
    /*Coarser neighbour - we need to find out the corresponding quadrant*/

    switch(dir) {
    case NTH:
      k = (i-1)/2; /*k is the corresponding quadrant of the neighbour in direction dir for cell i*/ 
      break;
    case STH:
      k = i/2; break;
    case EST:
      if(i <= 3)
	k = i-2;
      else k = i-4;
      break;
    case WST:
      if(i <= 1)
	k = i;
      else k = i-2;
      break;
    case LWR:
      k = i; break;
    case UPR:
      k = i-4; break;
    }
			      
    C -> face_neighbours[dir][0] -> flux_area[opp_dir][k] = 
      ((C -> flux_area[dir][0]) + (C -> flux_area[dir][1]) + (C -> flux_area[dir][2]) + (C -> flux_area[dir][3]));
	
#if AXI_OK	
    if((twoD == 2) && ((dir == EST) || (dir == WST)) && (k <= 1)) {
      if((C -> flux_area[dir][0] + C -> flux_area[dir][1]) > 0) {
	C -> face_neighbours[dir][0] -> face_r_centroid[opp_dir][k] = 
	  ((C -> face_r_centroid[dir][0])*(C -> flux_area[dir][0]) + (C -> face_r_centroid[dir][1])*(C -> flux_area[dir][1]))/
	  (C -> flux_area[dir][0] + C -> flux_area[dir][1]);
      } 
      else C -> face_neighbours[dir][0] -> face_r_centroid[opp_dir][k] = 0;      
    }
#endif
	  
    /*Flux area of quadrant equals sum of all of corresponding cell's quadrants (the face area)*/
  }
  else if(C -> face_neighbours[dir][0] -> cell_level == C -> cell_level) {
	    
    for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) 
      C -> face_neighbours[dir][0] -> flux_area[opp_dir][k] = C -> flux_area[dir][k];    

#if AXI_OK
    if((twoD == 2) && ((dir == EST) || (dir == WST))) {
      C -> face_neighbours[dir][0] -> face_r_centroid[opp_dir][0] = C -> face_r_centroid[dir][0];
      C -> face_neighbours[dir][0] -> face_r_centroid[opp_dir][1] = C -> face_r_centroid[dir][1];
    }
#endif
  }
  else if(C -> face_neighbours[dir][0] -> cell_level == (C -> cell_level + 1)) { 

    /*Can set neighbour's areas if the child's area either completely open or closed*/

    if(twoD != 0) {
      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) {
	if(
#if FLOAT_EQ_STABLE
	   EQ(C -> flux_area[dir][k], child_edge_length/2.0) || EQ(C -> flux_area[dir][k], 0)
#else
	   C -> flux_area[dir][k] == (child_edge_length/2.0) || (C -> flux_area[dir][k] == 0)
#endif
	   ) {
	  C -> face_neighbours[dir][k] -> flux_area[opp_dir][0] = (C -> flux_area[dir][k])/2.0;
	  C -> face_neighbours[dir][k] -> flux_area[opp_dir][1] = (C -> flux_area[dir][k])/2.0;
	  C -> face_neighbours[dir][k] -> flux_area[opp_dir][2] = 0;
	  C -> face_neighbours[dir][k] -> flux_area[opp_dir][3] = 0; /*Quadrants 2 and 3 should be 0*/

#if AXI_OK
	  if((twoD == 2) && ((dir == EST) || (dir == WST)) && (k < 2)) {
	    if(
#if FLOAT_EQ_STABLE
	       EQ(C -> flux_area[dir][k], 0)
#else
	       (C -> flux_area[dir][k] == 0)
#endif
	       ) { /*Fully obstructed (face_r_centroid would be 0)*/
	      C -> face_neighbours[dir][k] -> face_r_centroid[opp_dir][0] = 0;
	      C -> face_neighbours[dir][k] -> face_r_centroid[opp_dir][1] = 0;
	    }
	    else { /*Fully open*/
	      C -> face_neighbours[dir][k] -> face_r_centroid[opp_dir][0] = (C -> face_r_centroid[dir][k]) - (C -> cell_length)/8.0;
	      C -> face_neighbours[dir][k] -> face_r_centroid[opp_dir][1] = (C -> face_r_centroid[dir][k]) + (C -> cell_length)/8.0;
	    }
	  }
#endif	  
	}
      }
    }
    else {
      for(k = 0; k < MAX_NUM_NEIGHBOURS; k++) { /*Cycle through quadrants of the child*/
      					  
	if(
#if FLOAT_EQ_STABLE
	   EQ(C -> flux_area[dir][k], (child_face_area/MAX_NUM_NEIGHBOURS)) || EQ(C -> flux_area[dir][k], 0)
#else
	   (C -> flux_area[dir][k] == (child_face_area/MAX_NUM_NEIGHBOURS)) || (C -> flux_area[dir][k] == 0)
#endif
	   ) {	  
	  C -> face_neighbours[dir][k] -> flux_area[opp_dir][0] = (C -> flux_area[dir][k])/MAX_NUM_NEIGHBOURS;
	  C -> face_neighbours[dir][k] -> flux_area[opp_dir][1] = (C -> flux_area[dir][k])/MAX_NUM_NEIGHBOURS;
	  C -> face_neighbours[dir][k] -> flux_area[opp_dir][2] = (C -> flux_area[dir][k])/MAX_NUM_NEIGHBOURS;
	  C -> face_neighbours[dir][k] -> flux_area[opp_dir][3] = (C -> flux_area[dir][k])/MAX_NUM_NEIGHBOURS;
	}					  
      }
    }
  }
}

/*------------------------------------------------------------------*/

/**\brief Check if a given point is behind anyone of the spherical detonation waves (described in Deton) at time t */

short int check_if_in_deton_spheres(double time, double pt[])
{
  int i;
  double rd, rp;

  for(i = 0; i < Deton.num_deton_pts; i++) { /*Go through each detonation point*/
    rp = SQR(pt[0] - Deton.dpts[i][0]) + SQR(pt[1] - Deton.dpts[i][1]) + SQR(pt[2] - Deton.dpts[i][2]); /*Square of distance from point
													  to centre of detonation*/
    rd = SQR(time*(Deton.dcjvels[i])); /*Square of distance from point to detonation wave*/

    if(rp <= rd) /*Point is behind this detonation wave*/
      return(TRUE);
  }

  return(FALSE); /*Point isn't inside any region behind detonation wave*/
}

/*------------------------------------------------------------------*/

/**\brief Go through list of cells and see if any cells with un_det TRUE can now set it to FALSE*/

void comp_new_deton_cells(double time, List_leaf Head, List_leaf Tail)
{
  List_leaf L = Head;
  Cart_cell C;

  while(L != NULL) {
    C = L -> cell_loc;

    if(C -> un_det == TRUE) {
      Deton.deton_active = TRUE;
      if(check_if_in_deton_spheres(time, C -> centroid) == TRUE) {
	C -> un_det = FALSE;

#if 0
	if((C->parent->centroid[0]==7.81250000e-3) && (C->parent->centroid[1]==-7.81250000e-03) && (C->parent->centroid[2]==1.17187500e-01)
	   && (C->parent->cell_level==6)) {
	  printf("So at time %e L6 parent's child %hd passed, centroid %9.8e %9.8e %9.8e\n",time,C->child_num,C->centroid[0],C->centroid[1],C->centroid[2]);
	  printf("Radius is at %9.8e, C's radius is %9.8e\n",time*6717.4,SQRT(SQR(C->centroid[0]) + SQR(C->centroid[1]) + SQR(C->centroid[2]+0.25)));
	}
#endif
      }
    }

    if(L == Tail) 
      break;
    else L = L -> next;
  }
}

/*------------------------------------------------------------------*/
