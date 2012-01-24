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

/**\file Source file on constructing and traversing ADTs specific to cartesian subcells - geometric
 interrogation is a separate module*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ov_kernel.h"
#include "ov_adts.h"

#define NOTINCLUDE 0

extern short int twoD;
extern Point_tnode ADT2Dpoint;
extern Point_tnode ADT3Dpoint;
extern short int area_subcell_num;
extern short int volume_subcell_num;
extern double root_shift[3];
extern double root_scale;
extern double twoD_thickness;
extern double root_bbox[2][3];

/*------------------------------------------------------------------*/

/**\brief 
   Traverse ADT of points in 3D cell volume; count no. of points inside a solid and determine if only 1 actual body exists in 3D volume.
   Only determines if points are inside any body, not which.
*/

void traverse_volume_ADT(double shift[], double scale, int *num_subcells_inside, Point_tnode tnode, List_bbox *Intersected_bbox,
			 Body Solid, Body *temp_body, short int *is_unique_body, double bbox[])
{ 
  short int i;
  double global_point[3]; /*Transformed point in global 3D space*/

  for(i = 0; i < 3; i++)
    global_point[i] = (tnode -> point[i])*scale + shift[i];

  if(Solid == NULL) /*Must use the list as more than 1 body intersects 3D volume*/
    {
      List_bbox L = (*Intersected_bbox);

      while(L != NULL) /*Go thru each list node and test for bounding box intersections, then for actual intersections*/
	{	  
	  for(i = 0; i < 6; i++)
	    {
	      if(i < 3) 
		{
		  if(global_point[i] < (L -> tnode -> bounding_box[i])) /*Outside bounding box*/ 
		    break;
		}
	      else
		{
		  if(global_point[i%3] > (L -> tnode -> bounding_box[i])) /*Outside bounding box*/
		    break;
		}
	    }
	 	  
	  if(i == 6) /*Point within bounding box - test if actually inside the body*/
	    {
	      i = polyhio(global_point, L -> tnode -> Solid);

	      if(i != OUT_POLYH) /*If not outside, then treat as obstructed*/	      
		{
		  if(((*temp_body) == NULL) && ((*is_unique_body) == TRUE))
		    (*temp_body) = (L -> tnode -> Solid); 
		  else if(((*temp_body) != NULL) && ((*is_unique_body) == TRUE))
		    {
		      if((*temp_body) != (L -> tnode -> Solid)) /*Compare bodies*/
			(*is_unique_body) = FALSE; /*Bodies aren't the same; *is_unique_body FALSE*/
		    }

		  (*num_subcells_inside)++;
		  break;
		}
	    }
	  
	  L = L -> next;
	} 
      
      if(tnode -> left != NULL) /*Test if left region overlaps any bounding boxes*/
	{ 	  
	  L = (*Intersected_bbox);
	  
	  while(L != NULL)
	    {
	      if(((((tnode -> left -> upper_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) >
		  (L -> tnode -> bounding_box[(tnode -> level)%3])) &&
		 ((((tnode -> left -> lower_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) <
		  (L -> tnode -> bounding_box[(tnode -> level)%3 + 3])))
		{ 		  
		  /*Traverse down left link*/
		  traverse_volume_ADT(shift, scale, &(*num_subcells_inside), tnode -> left, &(*Intersected_bbox), Solid, 
				      &(*temp_body), &(*is_unique_body), bbox);
		  break;
		}	    	    
	      L = L -> next;
	    }
	}
      
      if(tnode -> right != NULL)
	{ /*Test if right region overlaps any bounding boxes*/
	  
	  L = (*Intersected_bbox);
	  
	  while(L != NULL)
	    {
	      if(((((tnode -> right -> upper_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) >
		  (L -> tnode -> bounding_box[(tnode -> level)%3])) &&
		 ((((tnode -> right -> lower_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) <
		  (L -> tnode -> bounding_box[(tnode -> level)%3 + 3])))
		{ 		  
		  /*Traverse down right link*/
		  traverse_volume_ADT(shift, scale, &(*num_subcells_inside), tnode -> right, &(*Intersected_bbox), Solid,
				      &(*temp_body), &(*is_unique_body), bbox);
		  break;
		}	    	    
	      L = L -> next;
	    }	  
	}  
    }
  else /*So need only use Solid pointer, not *Intersected_bbox*/
    {      
      if(*is_unique_body == -1) /*We're looking at IC conditions where a body ADT wasn't allocated - use bbox directly*/
	{
	  for(i = 0; i < 6; i++) 
	    {
	      if(i < 3)
		{
		  if(global_point[i] < bbox[i]) 
		    break;
		}
	      else
		{
		  if(global_point[i%3] > bbox[i])
		    break;
		}
	    }
	}
      else 
	{
	  for(i = 0; i < 6; i++)
	    {
	      if(i < 3) 
		{
		  if(global_point[i] < (Solid -> tnode -> bounding_box[i])) 
		    break;
		}
	      else
		{
		  if(global_point[i%3] > (Solid -> tnode -> bounding_box[i]))
		    break;
		}
	    }
	}
      
      if(i == 6) 
	{ 
	  i = polyhio(global_point, Solid); 

	  if(i != OUT_POLYH) 
	    (*num_subcells_inside)++; 
	}
      
      if(*is_unique_body == -1)
	{
	  if(tnode -> left != NULL)
	    { 	  
	      if(((((tnode -> left -> upper_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) >
		  bbox[(tnode -> level)%3]) &&
		 ((((tnode -> left -> lower_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) <
		  bbox[(tnode -> level)%3 + 3]))
		traverse_volume_ADT(shift, scale, &(*num_subcells_inside), tnode -> left, &(*Intersected_bbox), Solid, 
				    &(*temp_body), &(*is_unique_body), bbox); /*Traverse down left link*/
	    }
      
	  if(tnode -> right != NULL)
	    { 	  
	      if(((((tnode -> right -> upper_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) >
		  bbox[(tnode -> level)%3]) &&
		 ((((tnode -> right -> lower_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) <
		  bbox[(tnode -> level)%3 + 3]))
		traverse_volume_ADT(shift, scale, &(*num_subcells_inside), tnode -> right, &(*Intersected_bbox), Solid,
				    &(*temp_body), &(*is_unique_body), bbox); /*Traverse down right link*/
	    }
	}
      else
	{
	  if(tnode -> left != NULL)
	    { 	  
	      if(((((tnode -> left -> upper_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) >
		  (Solid -> tnode -> bounding_box[(tnode -> level)%3])) &&
		 ((((tnode -> left -> lower_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) <
		  (Solid -> tnode -> bounding_box[(tnode -> level)%3 + 3])))
		traverse_volume_ADT(shift, scale, &(*num_subcells_inside), tnode -> left, &(*Intersected_bbox), Solid, 
				    &(*temp_body), &(*is_unique_body), bbox); /*Traverse down left link*/
	    }
      
	  if(tnode -> right != NULL)
	    { 	  
	      if(((((tnode -> right -> upper_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) >
		  (Solid -> tnode -> bounding_box[(tnode -> level)%3])) &&
		 ((((tnode -> right -> lower_bound[(tnode -> level)%3]))*scale + shift[(tnode -> level)%3]) <
		  (Solid -> tnode -> bounding_box[(tnode -> level)%3 + 3])))
		traverse_volume_ADT(shift, scale, &(*num_subcells_inside), tnode -> right, &(*Intersected_bbox), Solid,
				    &(*temp_body), &(*is_unique_body), bbox); /*Traverse down right link*/
	    }
	}
    }
}

/*------------------------------------------------------------------*/

/**\brief 
   Traverse ADT of points in 2D plane - project 3D bounding boxes onto plane, then test if planar points inside bodies*/

void traverse_area_ADT(short int normal_direction, double shift[], double shift2D[], double plane_loc, double scale, 
		       int num_subcells_inside[], Point_tnode tnode, List_bbox *Intersected_bbox, Body Solid,
		       double r_solid[], int num2D_subcells_inside[])
{ 
  short int i;

  /*First transform point in the plane to a 3D point in global space*/

  double point3D[3]; /*Planar point in 3D space*/
  
  if(normal_direction == 1) /*North-south direction*/ 
    {
      point3D[0] = (tnode -> point[0])*scale + shift[0]; 
      point3D[1] = shift[1]; 
      point3D[2] = (tnode -> point[1])*scale + shift[2];

      /*point3D[1] is 'Y' co-ordinate of all; all points have same Y co-ordinate*/
    }
  else if(normal_direction == 2) /*Upper-lower direction*/
    {
      point3D[0] = (tnode -> point[0])*scale + shift[0]; 
      point3D[1] = (tnode -> point[1])*scale + shift[1];
      point3D[2] = shift[2];
    }  
  else /*East/West*/
    {
      point3D[0] = shift[0]; 
      point3D[1] = (tnode -> point[0])*scale + shift[1]; 
      point3D[2] = (tnode -> point[1])*scale + shift[2];
    }

  if(Solid == NULL) /*Use list*/
    {
      List_bbox L = (*Intersected_bbox);

      while(L != NULL) 
	{
	  for(i = 0; i < 6; i++)
	    {
	      if(i < 3)
		{
		  if(point3D[i] < (L -> tnode -> bounding_box[i]))
		    break; /*Point smaller than lower bound, so exit*/
		}
	      else
		{
		  if(point3D[i%3] > (L -> tnode -> bounding_box[i]))
		    break; /*Point larger than upper bound, so exit*/
		}
	    }

	  if(i == 6) 
	    {
	      i = polyhio(point3D, L -> tnode -> Solid);

	      if(i != OUT_POLYH)
		{
		  /*Quadrants numbered differently depending on direction*/
		  
		  if(normal_direction == 2) 
		    {
		      if((tnode -> point[0] < 0.5) && (tnode -> point[1] < 0.5)) { /*0th quadrant*/
			(num_subcells_inside[0])++;
#if AXI_OK
			if(twoD == 2)  
			  r_solid[0] += point3D[1];
#endif
		      }
		      else if((tnode -> point[0] < 0.5) && (tnode -> point[1] > 0.5)) { /*1st quadrant*/
			(num_subcells_inside[1])++;
#if AXI_OK
			if(twoD == 2)
			  r_solid[1] += point3D[1];
#endif
		      }
		      else if((tnode -> point[0] > 0.5) && (tnode -> point[1] < 0.5)) { /*2nd quadrant*/
			(num_subcells_inside[2])++;
#if AXI_OK
			if(twoD == 2)
			  r_solid[2] += point3D[1];
#endif
		      }
		      else {
			(num_subcells_inside[3])++; /*3rd quadrant*/
#if AXI_OK
			if(twoD == 2)
			  r_solid[3] += point3D[1];
#endif
		      }
		    }
		  else /*North, South, East, West faces*/
		    {
		      if((tnode -> point[0] < 0.5) && (tnode -> point[1] < 0.5)) { /*0th quadrant*/
			(num_subcells_inside[0])++;
#if AXI_OK
			if((twoD == 2) && (normal_direction == 0) && (point3D[2] < twoD_thickness)) {  /*r_solid only needed for east-west 
												     directions (quadrants 0 and 1 only)*/
			  r_solid[0] += point3D[1]; /*We want later to get average interface r co-ordinate*/
			  (num2D_subcells_inside[0])++;
			}
#endif
		      }
		      else if((tnode -> point[0] > 0.5) && (tnode -> point[1] < 0.5)) { /*1st quadrant*/
			(num_subcells_inside[1])++;
#if AXI_OK
			if((twoD == 2) && (normal_direction == 0) && (point3D[2] < twoD_thickness)) {
			  r_solid[1] += point3D[1];
			  (num2D_subcells_inside[1])++;
			}
#endif			
		      }
		      else if((tnode -> point[0] < 0.5) && (tnode -> point[1] > 0.5)) /*2nd quadrant*/
			(num_subcells_inside[2])++;
		      else (num_subcells_inside[3])++; /*3rd quadrant*/
		    }
		  
		  break;
		}
	    }

	  L = L -> next;
	}


      if(tnode -> left != NULL)
	{
	  L = (*Intersected_bbox);

	  while(L != NULL)
	    {
	      if((plane_loc >= (L -> tnode -> bounding_box[normal_direction])) && 
		 (plane_loc <= (L -> tnode -> bounding_box[normal_direction+3])))     
		{ /*Plane overlaps body in normal direction (flush overlaps permitted)*/
		  
		  if((((tnode -> left -> upper_bound[(tnode -> level)%2])*scale + shift2D[(tnode -> level)%2]) >
		      (L -> tnode -> plane_bbox[normal_direction][(tnode -> level)%2])) &&
		     (((tnode -> left -> lower_bound[(tnode -> level)%2])*scale + shift2D[(tnode -> level)%2]) <
		      (L -> tnode -> plane_bbox[normal_direction][(tnode -> level)%2 + 2])))		 
		    { /*Left region overlaps with projected bounding box*/

		      traverse_area_ADT(normal_direction, shift, shift2D, plane_loc, scale, num_subcells_inside, 
					tnode -> left, &(*Intersected_bbox), Solid, r_solid, num2D_subcells_inside); /*Traverse down left link*/
		      break;
		    } 
		}

	      L = L -> next;
	    }
	}

      if(tnode -> right != NULL)
	{
	  L = (*Intersected_bbox);

	  while(L != NULL)
	    {
	      if((plane_loc >= (L -> tnode -> bounding_box[normal_direction])) && 
		 (plane_loc <= (L -> tnode -> bounding_box[normal_direction+3])))     
		{ /*Overlapping in normal direction exists*/

		  if((((tnode -> right -> upper_bound[(tnode -> level)%2])*scale + shift2D[(tnode -> level)%2]) >
		      (L -> tnode -> plane_bbox[normal_direction][(tnode -> level)%2])) &&
		     (((tnode -> right -> lower_bound[(tnode -> level)%2])*scale + shift2D[(tnode -> level)%2]) <
		      (L -> tnode -> plane_bbox[normal_direction][(tnode -> level)%2 + 2])))		 
		    { /*Right region overlaps with projected bounding box*/

		      traverse_area_ADT(normal_direction, shift, shift2D, plane_loc, scale, num_subcells_inside, 
					tnode -> right, &(*Intersected_bbox), Solid, r_solid, num2D_subcells_inside); /*Traverse down right link*/

		      break;
		    } 
		}

	      L = L -> next;
	    }
	}
    }
  else /*Can just use Solid instead of *Intersected_bbox*/
    {
      for(i = 0; i < 6; i++)
	{
	  if(i < 3)
	    {
	      if(point3D[i] < (Solid -> tnode -> bounding_box[i]))
		break; /*Point smaller than lower bound, so exit*/
	    }
	  else
	    {
	      if(point3D[i%3] > (Solid -> tnode -> bounding_box[i]))
		break; /*Point larger than upper bound, so exit*/
	    }
	}
      
      if(i == 6) 
	{
	  i = polyhio(point3D, Solid);

	  if(i != OUT_POLYH)
	    {
	      if(normal_direction == 2) 
		{
		  if((tnode -> point[0] < 0.5) && (tnode -> point[1] < 0.5)) { /*0th quadrant*/
		    (num_subcells_inside[0])++;
#if AXI_OK		    
		    if(twoD == 2)
		      r_solid[0] += point3D[1];
#endif
		  }
		  else if((tnode -> point[0] < 0.5) && (tnode -> point[1] > 0.5)) { /*1st quadrant*/
		    (num_subcells_inside[1])++;
#if AXI_OK		    
		    if(twoD == 2)
		      r_solid[1] += point3D[1];
#endif
		  }
		  else if((tnode -> point[0] > 0.5) && (tnode -> point[1] < 0.5)) { /*2nd quadrant*/
		    (num_subcells_inside[2])++;
#if AXI_OK		    
		    if(twoD == 2)
		      r_solid[2] += point3D[1];
#endif
		  }
		  else {
		    (num_subcells_inside[3])++; /*3rd quadrant*/
#if AXI_OK		    
		    if(twoD == 2)
		      r_solid[3] += point3D[1];
#endif
		  }
		}
	      else /*North, South, East, West faces*/
		{
		  if((tnode -> point[0] < 0.5) && (tnode -> point[1] < 0.5)) { /*0th quadrant*/
		    (num_subcells_inside[0])++;
#if AXI_OK		    
		    if((twoD == 2) && (normal_direction == 0) && (point3D[2] < twoD_thickness)) {
		      r_solid[0] += point3D[1];
		      (num2D_subcells_inside[0])++;
		    }
#endif
		  }
		  else if((tnode -> point[0] > 0.5) && (tnode -> point[1] < 0.5)) { /*1st quadrant*/
		    (num_subcells_inside[1])++;
#if AXI_OK		    
		    if((twoD == 2) && (normal_direction == 0) && (point3D[2] < twoD_thickness)) {
		      r_solid[1] += point3D[1];
		      (num2D_subcells_inside[1])++;
		    }
#endif
		  }
		  else if((tnode -> point[0] < 0.5) && (tnode -> point[1] > 0.5)) /*2nd quadrant*/
		    (num_subcells_inside[2])++;
		  else (num_subcells_inside[3])++; /*3rd quadrant*/
		}
	    }
	}
      
      if(tnode -> left != NULL)
	{
	  if((plane_loc >= (Solid -> tnode -> bounding_box[normal_direction])) && 
	     (plane_loc <= (Solid -> tnode -> bounding_box[normal_direction+3])))     
	    { /*Overlapping in normal direction exists*/
			      	      
	      if((((tnode -> left -> upper_bound[(tnode -> level)%2])*scale + shift2D[(tnode -> level)%2]) >
		  (Solid -> tnode -> plane_bbox[normal_direction][(tnode -> level)%2])) &&
		 (((tnode -> left -> lower_bound[(tnode -> level)%2])*scale + shift2D[(tnode -> level)%2]) <
		  (Solid -> tnode -> plane_bbox[normal_direction][(tnode -> level)%2 + 2])))		 
		traverse_area_ADT(normal_direction, shift, shift2D, plane_loc, scale, 
				  num_subcells_inside, tnode -> left, &(*Intersected_bbox), Solid, r_solid, num2D_subcells_inside);
	    }
	}

      if(tnode -> right != NULL)
	{
	  if((plane_loc >= (Solid -> tnode -> bounding_box[normal_direction])) && 
	     (plane_loc <= (Solid -> tnode -> bounding_box[normal_direction+3])))     
	    { /*Overlapping in normal direction exists*/
	     	     	      	      	      
	      if((((tnode -> right -> upper_bound[(tnode -> level)%2])*scale + shift2D[(tnode -> level)%2]) >
		  (Solid -> tnode -> plane_bbox[normal_direction][(tnode -> level)%2])) &&
		 (((tnode -> right -> lower_bound[(tnode -> level)%2])*scale + shift2D[(tnode -> level)%2]) <
		  (Solid -> tnode -> plane_bbox[normal_direction][(tnode -> level)%2 + 2])))		 
		traverse_area_ADT(normal_direction, shift, shift2D, plane_loc, scale, 
				  num_subcells_inside, tnode -> right, &(*Intersected_bbox), Solid, r_solid, num2D_subcells_inside);
	    }	    
	}
    }
}

/*------------------------------------------------------------------*/

/**\brief Builds ADT of points in unit cube (corresp. to volume subcells)*/

short int build_volume_ADT(void)
{ 
  short int i, j, k; 
  int n = 0; /*No. subcell centroids scanned in*/
  double point[3]; 
  double opposite_point[3]; /*Diagonally opposite point (helps to get balanced ADT)*/

  double subcell_length = 1.0/volume_subcell_num;
  j = 0; k = 0;
  for(i = 0; i < volume_subcell_num/2; i++) 
    {
      point[0] = i*subcell_length + subcell_length/2.0;
      opposite_point[0] = (volume_subcell_num - 1 - i)*subcell_length + subcell_length/2.0; 

      for(j = 0; j < volume_subcell_num; j++) 
	{
	  point[1] = j*subcell_length + subcell_length/2.0; 
	  opposite_point[1] = (volume_subcell_num - 1 - j)*subcell_length + subcell_length/2.0;

	  for(k = 0; k < volume_subcell_num; k++) 
	    {
	      point[2] = k*subcell_length + subcell_length/2.0;
	      opposite_point[2] = (volume_subcell_num - 1 - k)*subcell_length + subcell_length/2.0;

	      n++;
	      
	      if(add_to_point_ADT(point, 3, &ADT3Dpoint, n) == ERROR) 
		break;
	      else 
		{
		  n++;

		  if(add_to_point_ADT(opposite_point, 3, &ADT3Dpoint, n) == ERROR) 
		    break;
		}
	    }
	}
    }

  if((i == volume_subcell_num/2) && (j == volume_subcell_num) && (k == volume_subcell_num))
    return(NOERROR); /*All subcells successfully put into ADT*/
  else return(ERROR);
}

/*------------------------------------------------------------------*/

/**\brief Builds ADT of points in unit square (corresp. to area subcells)*/

short int build_area_ADT(void)
{
  short int i, j; 
  int n = 0; /*No. of subcell centroids scanned in*/
  double point[2]; 
  double opposite_point[2]; /*Diagonally opposite point (helps to get balanced ADT)*/

  double subcell_length = 1.0/area_subcell_num;
  j = 0;
  for(i = 0; i < area_subcell_num/2; i++) 
    {
      point[0] = i*subcell_length + subcell_length/2.0; 
      opposite_point[0] = (area_subcell_num - 1 - i)*subcell_length + subcell_length/2.0; 

      for(j = 0; j < area_subcell_num; j++) 
	{
	  point[1] = j*subcell_length + subcell_length/2.0; 
	  opposite_point[1] = (area_subcell_num - 1 - j)*subcell_length + subcell_length/2.0; 

	  n++;
	  	  
	  if(add_to_point_ADT(point, 2, &ADT2Dpoint, n) == ERROR)
	    break;
	  else 
	    {
	      n++;
	      
	      if(add_to_point_ADT(opposite_point, 2, &ADT2Dpoint, n) == ERROR)
		break;
	    }
	}
    }
  
  if((i == area_subcell_num/2) && (j == area_subcell_num))
    return(NOERROR); /*All subcells successfully put into ADT*/
  else return(ERROR);
}

/*------------------------------------------------------------------*/

/**\brief 
   Adds point to ADT of points in plane/volume (global ADT2Dpoint/ADT3Dpoint).*/

short int add_to_point_ADT(double point[], short int dim, Point_tnode *tnode, int point_num)
{ 
  short int i;

  if((*tnode) == NULL) 
    {
      (*tnode) = malloc(sizeof(struct point_tnode));

      if((*tnode) != NULL)
	{ 
	  for(i = 0; i < dim; i++)
	    (*tnode) -> point[i] = point[i];
	  
	  /*Note for 2D points, (*tnode) -> point[2] is not set as in plane*/

	  (*tnode) -> level = 0;

	  for(i = 0; i < dim; i++)
	    {
	      (*tnode) -> lower_bound[i] = 0.0;
	      (*tnode) -> upper_bound[i] = 1.0; /*Bounds for unit hyperrectangle*/
	    }

	  (*tnode) -> left = NULL; (*tnode) -> right = NULL;

	  (*tnode) -> parent = NULL; (*tnode) -> point_num = point_num;

	  return(NOERROR); 
	}
      else return(ERROR);
    }
  else
    { /*Node already exists*/

      if(point[((*tnode) -> level)%dim] < 
	  0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%dim]) + ((*tnode) -> upper_bound[((*tnode) -> level)%dim])))
	{ /*Point lies in left subregion*/

	  if((*tnode) -> left == NULL)
	    {
	      (*tnode) -> left = malloc(sizeof(struct point_tnode));

	      if((*tnode) -> left != NULL)
		{		  
		  for(i = 0; i < dim; i++)
		    (*tnode) -> left -> point[i] = point[i];
		  
		  (*tnode) -> left -> level = (*tnode) -> level + 1;
		  
		  (*tnode) -> left -> left = NULL; (*tnode) -> left -> right = NULL;
		  
		  for(i = 0; i < dim; i++)
		    {
		      (*tnode) -> left -> lower_bound[i] = (*tnode) -> lower_bound[i];
		      (*tnode) -> left -> upper_bound[i] = (*tnode) -> upper_bound[i];
		    }

		  /*The upper bound has changed - now the midpoint of this subregion on this axis*/
		  
		  (*tnode) -> left -> upper_bound[((*tnode) -> level)%dim] = 
		    0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%dim]) + ((*tnode) -> upper_bound[((*tnode) -> level)%dim]));

		  (*tnode) -> left -> parent = (*tnode); (*tnode) -> left -> point_num = point_num;
		  
		  return(NOERROR); 
		} 
	      else return(ERROR);
	    }
	  else 
	    {
	      if(add_to_point_ADT(point, dim, &((*tnode) -> left), point_num) == ERROR)
		return(ERROR);
	      else return(NOERROR);
	    }	  
	}
      else
	{  /*Point lies in right subregion*/

	  if((*tnode) -> right == NULL)
	    {
	      (*tnode) -> right = malloc(sizeof(struct point_tnode));
	      
	      if((*tnode) -> right != NULL)
		{		  
		  for(i = 0; i < dim; i++)
		    (*tnode) -> right -> point[i] = point[i];
		  
		  (*tnode) -> right -> level = (*tnode) -> level + 1;
		  
		  (*tnode) -> right -> left = NULL; (*tnode) -> right -> right = NULL;

		  for(i = 0; i < dim; i++)
		    {
		      (*tnode) -> right -> lower_bound[i] = (*tnode) -> lower_bound[i];
		      (*tnode) -> right -> upper_bound[i] = (*tnode) -> upper_bound[i];
		    }
		  
		  /*The lower bound has changed - now the midpoint of this subregion on this axis*/
		  
		  (*tnode) -> right -> lower_bound[((*tnode) -> level)%dim] = 
		    0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%dim]) + ((*tnode) -> upper_bound[((*tnode) -> level)%dim]));

		  (*tnode) -> right -> parent = (*tnode); (*tnode) -> right -> point_num = point_num;
		  		  
		  return(NOERROR); 
		} 
	      else return(ERROR); 
	    }
	  else 
	    {
	      if(add_to_point_ADT(point, dim, &((*tnode) -> right), point_num) == ERROR)
		return(ERROR);
	      else return(NOERROR);
	    }
	}      
    }
}

/*------------------------------------------------------------------*/
