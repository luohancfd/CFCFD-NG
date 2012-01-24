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



/**\file Source file for constructing and using object ADTs for point inclusion interrogation.  Please see
 gio_kernel.c for source on point inclusion queries*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gio_kernel.h"
#include "gio_adts.h"
#include "gio_lists.h"

#define NOTINCLUDE 0

extern Body_tnode ADTbody;
extern double root_3D_bbox[2][3];

/*------------------------------------------------------------------*/

/**\brief Build ADT of bodies (just rectangular prisms for the moment).  
   Ignores bounding boxes outside domain, trims all bounding boxes so they lie within root.
   Also calculates projected 6D bounding boxes in 4D hyperspace.  Returns no. bodies added to ADT, else ERROR.
*/

int build_body_ADT() 
{                             
  FILE *fin;
  char buffer[512];
  short int i;
  short int stop_scanning = FALSE;
  
  int n = 0; /*No. bodies scanned*/

  double lb[3]; /*Lower bound/corner of bounding box*/
  double ub[3]; /*Upper bound/corner of bounding box*/
  
  double bbox[6]; /*Bounding box in 6D*/
  double plane_bbox[3][4]; /*4D Projected bounding box on plane*/
    
  double root_lb[6]; double root_ub[6]; 

  Body Somebody;

  for(i = 0; i < 6; i++)
    {
      root_lb[i] = root_3D_bbox[0][i%3];
      root_ub[i] = root_3D_bbox[1][i%3];
    }
  
  fin = fopen("ov_bodies.par","r"); /*Standard input file for storing bodies*/

  if(fin == NULL)
    {
      printf("build_body_ADT() could not open file ov_bodies.par\n");
      return(ERROR);
    }

  /*Ignore first 3 lines*/
  fgets(buffer, sizeof(buffer), fin);
  fgets(buffer, sizeof(buffer), fin);
  fgets(buffer, sizeof(buffer), fin);

  while(stop_scanning == FALSE)
    {
      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "lb: %lf", &lb[0]) == 0) /*No more bodies to read*/
	{
	  printf("build_body_ADT() has no more correct bodies to read\n");
	  stop_scanning = TRUE;
	}
      else if(sscanf(buffer, "lb: %lf %lf %lf", &lb[0], &lb[1], &lb[2]) < 3)
	{
	  printf("build_body_ADT() read lower bounds incorrectly\n");
	  return(ERROR);
	}
      else /*First line read correctly*/
	{
	  fgets(buffer, sizeof(buffer), fin);

	  if(sscanf(buffer, "ub: %lf %lf %lf", &ub[0], &ub[1], &ub[2]) < 3)
	    {
	      printf("build_body_ADT() read upper bounds incorrectly\n");
	      return(ERROR);
	    }
	  else if((ub[0] < lb[0]) || (ub[1] < lb[1]) || (ub[2] < lb[2]))
	    {
	      printf("build_body_ADT() saw a weird body for body_num %d.  Please check it.\n", n+1);
	      return(ERROR);
	    } 
	  else /*Body has been read in correctly*/
	    {
	      fgets(buffer, sizeof(buffer), fin); /*Space in between bodies*/
		      
	      for(i = 0; i < 6; i++)
		{
		  if(i < 3)
		    bbox[i] = lb[i];
		  else
		    bbox[i] = ub[i%3];
		}

	      if(((bbox[0] <= root_3D_bbox[1][0]) && (bbox[1] <= root_3D_bbox[1][1]) && (bbox[2] <= root_3D_bbox[1][2])) &&
		 ((bbox[3] >= root_3D_bbox[0][0]) && (bbox[4] >= root_3D_bbox[0][1]) && (bbox[5] >= root_3D_bbox[0][2])))
		{ 
		  /*Bounding box intersects octree root's bounding box*/

		  for(i = 0; i < 6; i++) /*Trim bounding boxes that exceed octree root's domain to domain limits*/
		    {			
		      if(i < 3)
			{
			  if(bbox[i] < root_3D_bbox[0][i])
			    bbox[i] = root_3D_bbox[0][i];
			}
		      else
			{
			  if(bbox[i] > root_3D_bbox[1][i%3])
			    bbox[i] = root_3D_bbox[1][i%3];
			}
		    }

		  /*Bounding box on x plane - is YZ plane, extract Y and Z co-ords of 3D bounding box on plane*/
		  plane_bbox[0][0] = bbox[1];
		  plane_bbox[0][1] = bbox[2];
		  plane_bbox[0][2] = bbox[4];
		  plane_bbox[0][3] = bbox[5];

		  /*Bounding box on y plane*/
		  plane_bbox[1][0] = bbox[0];
		  plane_bbox[1][1] = bbox[2];
		  plane_bbox[1][2] = bbox[3];
		  plane_bbox[1][3] = bbox[5];

		  /*Bounding box on z plane*/
		  plane_bbox[2][0] = bbox[0];
		  plane_bbox[2][1] = bbox[1];
		  plane_bbox[2][2] = bbox[3];
		  plane_bbox[2][3] = bbox[4];

		  Somebody = malloc(sizeof(struct body));
		  		  
		  if(add_to_body_ADT(Somebody, bbox, plane_bbox, &ADTbody, n+1) == NOERROR)
		    n++;

		  Somebody = NULL;
		} 
	    }
	}
    }

  fclose(fin);
  
  return(n);    
}

/*------------------------------------------------------------------*/

/**\brief Build ADT of faces given the polygonal facets of a polyhedron body*/

short int build_face_ADT(Body Solid, short int compute_xy_bounds)
{
  Body2D Face;
  int i, j;
  double bbox[4];
  double root_bounds[2][2]; 

  Solid -> Faces = NULL; /*Safely initialize*/

  /*If required, get minimum and maximum bounds (on XY plane) for solid - represent as point in 4D space*/

  if(compute_xy_bounds != FALSE) 
    {
      root_bounds[0][0] = Solid -> points[0][0];
      root_bounds[0][1] = Solid -> points[0][1]; /*First row is lower bound*/
      root_bounds[1][0] = root_bounds[0][0];
      root_bounds[1][1] = root_bounds[0][1]; /*Second row is upper bound*/

      for(i = 0; i < Solid -> num_points; i++)
	{
	  if(Solid -> points[i][0] < root_bounds[0][0])
	    root_bounds[0][0] = Solid -> points[i][0];
	  else if(Solid -> points[i][0] > root_bounds[1][0])
	    root_bounds[1][0] = Solid -> points[i][0]; /*X co-ordinate*/

	  if(Solid -> points[i][1] < root_bounds[0][1])
	    root_bounds[0][1] = Solid -> points[i][1];
	  else if(Solid -> points[i][1] > root_bounds[1][1])
	    root_bounds[1][1] = Solid -> points[i][1]; /*Y co-ordinate*/
	}

      Solid -> xy_bound[0] = root_bounds[0][0];
      Solid -> xy_bound[1] = root_bounds[0][1];
      Solid -> xy_bound[2] = root_bounds[1][0];
      Solid -> xy_bound[3] = root_bounds[1][1];  
    }
  else 
    {
      root_bounds[0][0] = Solid -> xy_bound[0];
      root_bounds[0][1] = Solid -> xy_bound[1];
      root_bounds[1][0] = Solid -> xy_bound[2];
      root_bounds[1][1] = Solid -> xy_bound[3];
    }

  /*Add faces to ADT of edges*/

  for(i = 0; i < Solid -> num_faces; i++) /*Go through all faces*/
    {
      Face = malloc(sizeof(struct body2D));
      if(Face == NULL)
	return(ERROR);

      Face -> points = malloc(sizeof(double *)*(Solid->face_size[i])); /*Define points array for this face*/
      if(Face -> points == NULL)
	return(ERROR);

      Face -> edges = malloc(sizeof(int *)*(Solid->face_size[i])); /*As many edges as there are points in a closed polygon*/
      if(Face -> edges == NULL)
	return(ERROR);
      
      Face -> num_points = Solid -> face_size[i];

      bbox[0] = Solid -> points[Solid->faces[i][0]][0];
      bbox[1] = Solid -> points[Solid->faces[i][0]][1];
      bbox[2] = bbox[0];
      bbox[3] = bbox[1];

      for(j = 0; j < (Solid -> face_size[i]); j++) /*Go through each point (and edge) on the face*/
	{
	  /*Again compute bounding box for this particular polygon*/

	  if(Solid -> points[Solid->faces[i][j]][0] < bbox[0])
	    bbox[0] = Solid -> points[Solid->faces[i][j]][0];
	  else if(Solid -> points[Solid->faces[i][j]][0] > bbox[2])
	    bbox[2] = Solid -> points[Solid->faces[i][j]][0];

	  if(Solid -> points[Solid->faces[i][j]][1] < bbox[1])
	    bbox[1] = Solid -> points[Solid->faces[i][j]][1];
	  else if(Solid -> points[Solid->faces[i][j]][1] > bbox[3])
	    bbox[3] = Solid -> points[Solid->faces[i][j]][1];

	  Face -> points[j] = malloc(sizeof(double)*3);
	  if(Face -> points[j] == NULL)
	    return(ERROR);

	  Face -> edges[j] = malloc(sizeof(int)*2); /*Edges only need start and end point*/
	  if(Face -> edges[j] == NULL)
	    return(ERROR);

	  /*Define points in order for the face*/

	  Face -> points[j][0] = Solid -> points[Solid->faces[i][j]][0];
	  Face -> points[j][1] = Solid -> points[Solid->faces[i][j]][1];
	  Face -> points[j][2] = Solid -> points[Solid->faces[i][j]][2];	  

	  /*Edges simply connect the points*/

	  if(j == (Solid -> face_size[i] - 1))
	    Face -> edges[j][1] = 0; /*End point of last edge is 0th point*/
	  else Face -> edges[j][1] = j+1;

	  Face -> edges[j][0] = j;
	}

      /*Points and edges on face now defined - build the edge ADT for the face*/

      if(build_edge_ADT(Face) == ERROR)
	{
	  printf("Build_edge_ADT(): Error.  Exiting ...\n");
	  return(ERROR);
	}

      /*Now can add this face to ADT of faces*/

      if(add_to_face_ADT(Face, bbox, root_bounds, &(Solid -> Faces), i) == ERROR)
	{
	  printf("Add_to_face_ADT(): Error.  Exiting ...\n");
	  return(ERROR);
	}
      
      Face = NULL;
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Build ADT of edges given a 2D polygon*/

short int build_edge_ADT(Body2D Poly)
{
  Edge3D Edge;
  int i;
  double bbox[2];
  double root_bounds[2][1];

  Poly -> Edges = NULL; /*Must safely initialize*/

  /*Get minimum and maximum bounds (on x-axis) for polygon*/

  root_bounds[0][0] = Poly -> points[0][0];
  root_bounds[1][0] = root_bounds[0][0];

  for(i = 0; i < Poly -> num_points; i++)
    {
      if(Poly -> points[i][0] < root_bounds[0][0])
	root_bounds[0][0] = Poly -> points[i][0];
      else if(Poly -> points[i][0] > root_bounds[1][0])
	root_bounds[1][0] = Poly -> points[i][0];
    }

  Poly -> x_bound[0] = root_bounds[0][0]; 
  Poly -> x_bound[1] = root_bounds[1][0];

  /*Add edges to ADT of edges (there are as many points as edges for non-degenerate polygons)*/

  for(i = 0; i < Poly -> num_points; i++)
    {
      Edge = malloc(sizeof(struct edge3D));
      if(Edge == NULL)
	return(ERROR);

      Edge -> Poly = Poly;

      Edge -> start[0] = Poly -> points[Poly -> edges[i][0]][0]; /*Edge matrix only references points which constitute its start/end*/
      Edge -> start[1] = Poly -> points[Poly -> edges[i][0]][1];
      Edge -> start[2] = Poly -> points[Poly -> edges[i][0]][2];
      Edge -> end[0] = Poly -> points[Poly -> edges[i][1]][0];
      Edge -> end[1] = Poly -> points[Poly -> edges[i][1]][1];
      Edge -> end[2] = Poly -> points[Poly -> edges[i][1]][2];

      bbox[0] = MIN(Edge -> start[0],Edge -> end[0]); /*Even projected onto x-axis, edges might be defined clockwise/counter-clockwise*/
      bbox[1] = MAX(Edge -> start[0],Edge -> end[0]);
      
      if(add_to_edge_ADT(Edge, bbox, root_bounds, &(Poly -> Edges), i) == ERROR)
	{
	  printf("Add_to_edge_ADT(): Error.  Exiting ...\n");
	  return(ERROR);
	}

      Edge = NULL;
    }

  /*Finally precompute partial scalar triple product for determining later a queried point's status relative to polygon plane*/

  if(precompute_normal(Poly) == ERROR) 
    {
      printf("Build_edge_ADT(): Polygon appears to be degenerate.  Exiting ...\n");
      return(ERROR);
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Traverse ADT of bodies' bounding boxes, return list of bounding boxes which intersect specified region*/

short int traverse_body_ADT(double lower_bound[], double upper_bound[], Body_tnode tnode, List_bbox *L) 
{ 
  if(tnode != NULL) /*If there aren't any bodies, there's no point in going further*/
    {
      short int i;

      for(i = 0; i < 6; i++)
	{
	  if(i < 3)
	    {
	      if((tnode -> bounding_box[i]) > upper_bound[i]) /*Body is outside region*/
		break;
	    }
	  else if((tnode -> bounding_box[i]) < lower_bound[i])
	    break;
	}
      
      if(i == 6) /*i == 6, then the body does intersect the region - we can add it to the list*/
	if(add_to_list_bbox(tnode, &(*L)) == ERROR)
	  return(ERROR);
            
      /*Now test the left and right links for intersections*/
		  
      if(tnode -> left != NULL)
	{
	  if((tnode -> left -> upper_bound[(tnode -> level)%6] >= lower_bound[(tnode -> level)%6]) &&
	     (tnode -> left -> lower_bound[(tnode -> level)%6] <= upper_bound[(tnode -> level)%6]))
	    if(traverse_body_ADT(lower_bound, upper_bound, tnode -> left, &(*L)) == ERROR)
	      return(ERROR);

	  /*Overlap in left subregion exists; traverse down left link*/	      
	}
	  
      if(tnode -> right != NULL)
	{
	  if((tnode -> right -> upper_bound[(tnode -> level)%6] >= lower_bound[(tnode -> level)%6]) &&
	     (tnode -> right -> lower_bound[(tnode -> level)%6] <= upper_bound[(tnode -> level)%6]))
	    if(traverse_body_ADT(lower_bound, upper_bound, tnode -> right, &(*L)) == ERROR)
	      return(ERROR);
	      
	  /*Overlap in right subregion exists; traverse down right link*/
	}	
    } 
  
  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Traverse ADT of edges in a polygon; for all intersection candidates calculate and sum global linking numbers*/

short int traverse_edge_ADT(double point[], double lower_bound[], double upper_bound[], Edge1D_tnode tnode, int *link_num, 
			    double *vtx_angle)
{
  if(tnode != NULL)
    {
      short int i;
            
      for(i = 0; i < 2; i++)
	{
	  if(i == 0)
	    {
	      if((tnode -> bounding_box[i]) > upper_bound[i]) 
		break;
	    }
	  else if((tnode -> bounding_box[i]) < lower_bound[i])
	    break;
	}

      if(i == 2) 
	{
	  i = get_link_num(tnode -> Edge, point, link_num);

	  if(i == ON_POLY_EDGE)
	    return(i);
	  else if(i == ON_POLY_VTX)
	    {
	      if((*vtx_angle) < 0) /*So let's get the absolute value of the inner angle of the projected polygon at this vertex*/
		{
		  if(point[0] == (tnode -> Edge -> start[0])) /*Point is at start of edge (enough to compare x co-ordinates)*/
		    { /*We want the previous edge - all polygons on the polyhedron will hopefully be 'well defined'*/
		      
		      if(tnode -> edge_num == 0) /*Is the first edge*/
			*vtx_angle = get_vtx_angle(tnode -> Edge -> Poly, (tnode -> Edge -> Poly -> num_points) - 1, 0);
		      else *vtx_angle = get_vtx_angle(tnode -> Edge -> Poly, (tnode -> edge_num) - 1, tnode -> edge_num);
		    }
		  else /*Point is at end of edge (the only option left)*/
		    { /*We want the next edge*/

		      if(tnode -> edge_num == (tnode -> Edge -> Poly -> num_points) - 1) /*Is the last edge*/
			*vtx_angle = get_vtx_angle(tnode -> Edge -> Poly, tnode -> edge_num, 0);
		      else *vtx_angle = get_vtx_angle(tnode -> Edge -> Poly, tnode -> edge_num, (tnode -> edge_num) + 1);
		    }
		}

	      return(i);
	    }
	}
      
      if(tnode -> left != NULL)
	{
	  if((tnode -> left -> upper_bound[(tnode -> level)%2] >= lower_bound[(tnode -> level)%2]) &&
	     (tnode -> left -> lower_bound[(tnode -> level)%2] <= upper_bound[(tnode -> level)%2]))
	    {
	      i = traverse_edge_ADT(point, lower_bound, upper_bound, tnode -> left, link_num, vtx_angle); 
	      
	      if((i == ON_POLY_VTX) || (i == ON_POLY_EDGE))
		return(i);
	    }
	}

      if(tnode -> right != NULL)
	{
	  if((tnode -> right -> upper_bound[(tnode -> level)%2] >= lower_bound[(tnode -> level)%2]) &&
	     (tnode -> right -> lower_bound[(tnode -> level)%2] <= upper_bound[(tnode -> level)%2]))
	    {
	      i = traverse_edge_ADT(point, lower_bound, upper_bound, tnode -> right, link_num, vtx_angle); 
	      
	      if((i == ON_POLY_VTX) || (i == ON_POLY_EDGE))
		return(i);
	    }
	}
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Traverse ADT of faces (projected onto XY plane) of polyhedron*/

short int traverse_face_ADT(double point[], double lower_bound[], double upper_bound[], Body2D_tnode tnode, double *link_num)
{
  if(tnode != NULL)
    {
      short int i;
      double j;
      double s;
      double vtx_angle = -1; /*We do want to know what this inner vertex angle is*/
      
      for(i = 0; i < 4; i++)
	{
	  if(i < 2)
	    {
	      if((tnode -> bounding_box[i]) > upper_bound[i]) /*Body is outside region*/
		break;
	    }
	  else if((tnode -> bounding_box[i]) < lower_bound[i])
	    break;
	}

      if(i == 4) /*Point potentially intersects this polygon - get linking number for this polygon*/ 
	{ 
	  i = polyio(point, tnode -> Poly, &vtx_angle);
	  
	  if(i != OUT_POLY)
	    {
	      if(i == IN_POLY) 
		s = 0.5;
	      else if(i == ON_POLY_EDGE)
		s = 0.25;
	      else s = (0.5*vtx_angle)/(2.0*PI); /*Projected onto a vertex*/
	      
	      j = compute_side_loc(point, tnode -> Poly); /*See if point is upside/downside of facet*/
	      
	      if((j == 0) && (tnode -> Poly -> norm[2] != 0))  /*Point is on a surface facet - ignore projections
								 which form semi-infinite planes parallel along z-axis*/
		return(ON_POLYH);
	      else (*link_num) += s*j;		
	    }
	}
   	
      if(tnode -> left != NULL)
	{
	  if((tnode -> left -> upper_bound[(tnode -> level)%4] >= lower_bound[(tnode -> level)%4]) &&
	     (tnode -> left -> lower_bound[(tnode -> level)%4] <= upper_bound[(tnode -> level)%4]))
	    {
	      i = traverse_face_ADT(point, lower_bound, upper_bound, tnode -> left, link_num);
	      if(i == ON_POLYH)
		return(i);
	    }
	}

      if(tnode -> right != NULL)
	{
	  if((tnode -> right -> upper_bound[(tnode -> level)%4] >= lower_bound[(tnode -> level)%4]) &&
	     (tnode -> right -> lower_bound[(tnode -> level)%4] <= upper_bound[(tnode -> level)%4]))
	    {
	      i = traverse_face_ADT(point, lower_bound, upper_bound, tnode -> right, link_num);
	      if(i == ON_POLYH)
		return(i);
	    }
	}
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Adds body (and its bounding box - represented as point in hypercube) to ADT of bodies (global ADTbody)*/   

short int add_to_body_ADT(Body Somebody, double bounding_box[], double plane_bbox[][4], Body_tnode *tnode, int body_num)
{ 
  short int i, j;

  if((*tnode) == NULL) 
    {
      (*tnode) = malloc(sizeof(struct body_tnode));
      
      if((*tnode) != NULL)
	{
	  (*tnode) -> Solid = Somebody;

	  (*tnode) -> Solid -> tnode = (*tnode); /*Set parameters in (*tnode) -> Solid*/

	  /*Now set variables in *tnode*/

	  (*tnode) -> level = 0;
	      
	  for(i = 0; i < 6; i++)
	    (*tnode) -> bounding_box[i] = bounding_box[i];

	  for(i = 0; i < 3; i++)
	    {
	      for(j = 0; j < 4; j++)
		(*tnode) -> plane_bbox[i][j] = plane_bbox[i][j];
	    }

	  /*Compute bounding box of whole region*/

	  for(i = 0; i < 6; i++)
	    {
	      (*tnode) -> lower_bound[i] = root_3D_bbox[0][i%3];
	      (*tnode) -> upper_bound[i] = root_3D_bbox[1][i%3];
	    }

	  (*tnode) -> left = NULL; (*tnode) -> right = NULL;

	  (*tnode) -> parent = NULL; (*tnode) -> body_num = body_num;

	  return(NOERROR); 
	}
      else return(ERROR); 
    }
  else /*Node already exists*/
    {      
      if(bounding_box[((*tnode) -> level)%6] < 
	  0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%6]) + ((*tnode) -> upper_bound[((*tnode) -> level)%6])))
	{ /*Overlap in left subregion of (*tnode -> level)%6 th axis*/
	  
	  if((*tnode) -> left == NULL) 
	    {
	      (*tnode) -> left = malloc(sizeof(struct body_tnode));

	      if((*tnode) -> left != NULL) 
		{
		  (*tnode) -> left -> Solid = Somebody;
	
		  (*tnode) -> left -> Solid -> tnode = (*tnode) -> left;
		     		      
		  (*tnode) -> left -> level = (*tnode) -> level + 1;
		      
		  for(i = 0; i < 6; i++)
		    {
		      (*tnode) -> left -> bounding_box[i] = bounding_box[i];
		      (*tnode) -> left -> lower_bound[i] = (*tnode) -> lower_bound[i]; /*Lower bound in hypercube (temporary value)*/
		      (*tnode) -> left -> upper_bound[i] = (*tnode) -> upper_bound[i]; /*Upper bound in hypercube*/
		    }

		  for(i = 0; i < 3; i++)
		    {
		      for(j = 0; j < 4; j++)
			(*tnode) -> left -> plane_bbox[i][j] = plane_bbox[i][j];
		    }

		  /*The upper bound has changed - now the midpoint of this subregion on this axis*/
		      
		  (*tnode) -> left -> upper_bound[((*tnode) -> level)%6] = 
		    0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%6]) + ((*tnode) -> upper_bound[((*tnode) -> level)%6]));     
		      		      
		  (*tnode) -> left -> left = NULL; (*tnode) -> left -> right = NULL;

		  (*tnode) -> left -> parent = (*tnode); (*tnode) -> left -> body_num = body_num;

		  return(NOERROR);
		}
	      else return(ERROR);
	    }
	  else /*Node already exists*/
	    {
	      if(add_to_body_ADT(Somebody, bounding_box, plane_bbox, &((*tnode) -> left), body_num) == ERROR)
		return(ERROR);
	      else return(NOERROR);
	    }	  
	}
      else /*Overlap in right subregion of (*tnode -> level)%6 th axis*/
	{ 
	  if((*tnode) -> right == NULL) 
	    {
	      (*tnode) -> right = malloc(sizeof(struct body_tnode));

	      if((*tnode) -> right != NULL)
		{
		  (*tnode) -> right -> Solid = Somebody;
		      
		  (*tnode) -> right -> Solid -> tnode = (*tnode) -> right;
		     
		  (*tnode) -> right -> level = (*tnode) -> level + 1;
		      		      	      
		  for(i = 0; i < 6; i++)
		    {
		      (*tnode) -> right -> bounding_box[i] = bounding_box[i];
		      (*tnode) -> right -> lower_bound[i] = (*tnode) -> lower_bound[i]; /*Lower bound in hypercube*/
		      (*tnode) -> right -> upper_bound[i] = (*tnode) -> upper_bound[i]; /*Upper bound in hypercube*/
		    }

		  for(i = 0; i < 3; i++)
		    {
		      for(j = 0; j < 4; j++)
			(*tnode) -> right -> plane_bbox[i][j] = plane_bbox[i][j];
		    }
		      
		  /*The lower bound has changed - now the midpoint of this subregion on this axis*/

		  (*tnode) -> right -> lower_bound[((*tnode) -> level)%6] = 
		    0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%6]) + ((*tnode) -> upper_bound[((*tnode) -> level)%6]));

		  (*tnode) -> right -> left = NULL; (*tnode) -> right -> right = NULL;

		  (*tnode) -> right -> parent = (*tnode); (*tnode) -> right -> body_num = body_num;

		  return(NOERROR);		    
		}
	      else return(ERROR);
	    }
	  else /*Node already exists*/
	    {
	      if(add_to_body_ADT(Somebody, bounding_box, plane_bbox, &((*tnode) -> right), body_num) == ERROR)
		return(ERROR);
	      else return(NOERROR);
	    }	  
	}      
    }
}

/*------------------------------------------------------------------*/

/**\brief Add node to ADT of projected edges (onto x-axis, for simplicity)*/

short int add_to_edge_ADT(Edge3D Someedge, double bounding_box[], double root_1D_bbox[][1], Edge1D_tnode *tnode, int edge_num)
{
  short int i;

  if((*tnode) == NULL)
    {
      (*tnode) = malloc(sizeof(struct edge1D_tnode));

      if((*tnode) != NULL)
	{
	  (*tnode) -> Edge = Someedge;

	  (*tnode) -> Edge -> tnode = (*tnode);

	  (*tnode) -> level = 0;

	  for(i = 0; i < 2; i++)
	    {
	      (*tnode) -> bounding_box[i] = bounding_box[i];
	      (*tnode) -> lower_bound[i] = root_1D_bbox[0][i%1];
	      (*tnode) -> upper_bound[i] = root_1D_bbox[1][i%1];
	    }

	  (*tnode) -> left = NULL; (*tnode) -> right = NULL;

	  (*tnode) -> parent = NULL; (*tnode) -> edge_num = edge_num;

	  return(NOERROR);
	}
      else return(ERROR);
    }
  else
    {
      if(bounding_box[((*tnode) -> level)%2] < 
	 0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%2]) + ((*tnode) -> upper_bound[((*tnode) -> level)%2])))
	{
	  if((*tnode) -> left == NULL)
	    {
	      (*tnode) -> left = malloc(sizeof(struct edge1D_tnode));

	      if((*tnode) -> left != NULL)
		{
		  (*tnode) -> left -> Edge = Someedge;

		  (*tnode) -> left -> Edge -> tnode = (*tnode) -> left;

		  (*tnode) -> left -> level = (*tnode) -> level + 1;

		  for(i = 0; i < 2; i++)
		    {
		      (*tnode) -> left -> bounding_box[i] = bounding_box[i];
		      (*tnode) -> left -> lower_bound[i] = (*tnode) -> lower_bound[i]; 
		      (*tnode) -> left -> upper_bound[i] = (*tnode) -> upper_bound[i]; 
		    }

		  (*tnode) -> left -> upper_bound[((*tnode) -> level)%2] = 
		    0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%2]) + ((*tnode) -> upper_bound[((*tnode) -> level)%2]));

		  (*tnode) -> left -> left = NULL; (*tnode) -> left -> right = NULL;

		  (*tnode) -> left -> parent = (*tnode); (*tnode) -> left -> edge_num = edge_num;

		  return(NOERROR);
		}
	      else return(ERROR);
	    }
	  else
	    {
	      if(add_to_edge_ADT(Someedge, bounding_box, root_1D_bbox, &((*tnode) -> left), edge_num) == ERROR)
		return(ERROR);
	      else return(NOERROR);
	    }
	}
      else
	{
	  if((*tnode) -> right == NULL) 
	    {
	      (*tnode) -> right = malloc(sizeof(struct edge1D_tnode));

	      if((*tnode) -> right != NULL)
		{
		  (*tnode) -> right -> Edge = Someedge;
		      
		  (*tnode) -> right -> Edge -> tnode = (*tnode) -> right;
		     
		  (*tnode) -> right -> level = (*tnode) -> level + 1;
		      		      	      
		  for(i = 0; i < 2; i++)
		    {
		      (*tnode) -> right -> bounding_box[i] = bounding_box[i];
		      (*tnode) -> right -> lower_bound[i] = (*tnode) -> lower_bound[i]; 
		      (*tnode) -> right -> upper_bound[i] = (*tnode) -> upper_bound[i]; 
		    }
		      
		  (*tnode) -> right -> lower_bound[((*tnode) -> level)%2] = 
		    0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%2]) + ((*tnode) -> upper_bound[((*tnode) -> level)%2]));

		  (*tnode) -> right -> left = NULL; (*tnode) -> right -> right = NULL;

		  (*tnode) -> right -> parent = (*tnode); (*tnode) -> right -> edge_num = edge_num;

		  return(NOERROR);		    
		}
	      else return(ERROR);
	    }
	  else /*Node already exists*/
	    {
	      if(add_to_edge_ADT(Someedge, bounding_box, root_1D_bbox, &((*tnode) -> right), edge_num) == ERROR)
		return(ERROR);
	      else return(NOERROR);
	    }	 
	}
    }
}

/*------------------------------------------------------------------*/

/**\brief Adds face of polyhedron to ADT of faces for polyhedron*/

short int add_to_face_ADT(Body2D Someface, double bounding_box[], double root_2D_bbox[][2], Body2D_tnode *tnode, int poly_num)
{
  short int i;

  if((*tnode) == NULL)
    {
      (*tnode) = malloc(sizeof(struct body2D_tnode));

      if((*tnode) != NULL)
	{
	  (*tnode) -> Poly = Someface;

	  (*tnode) -> Poly -> tnode = (*tnode);

	  (*tnode) -> level = 0;

	  for(i = 0; i < 4; i++)
	    { 
	      (*tnode) -> bounding_box[i] = bounding_box[i];
	      (*tnode) -> lower_bound[i] = root_2D_bbox[0][i%2];
	      (*tnode) -> upper_bound[i] = root_2D_bbox[1][i%2];
	    }

	  (*tnode) -> left = NULL; (*tnode) -> right = NULL;

	  (*tnode) -> parent = NULL; (*tnode) -> poly_num = poly_num;

	  return(NOERROR);
	}
      else return(ERROR);
    }
  else
    {
      if(bounding_box[((*tnode) -> level)%4] < 
	 0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%4]) + ((*tnode) -> upper_bound[((*tnode) -> level)%4])))
	{ 	  
	  if((*tnode) -> left == NULL) 
	    {
	      (*tnode) -> left = malloc(sizeof(struct body2D_tnode));

	      if((*tnode) -> left != NULL) 
		{
		  (*tnode) -> left -> Poly = Someface;
	
		  (*tnode) -> left -> Poly -> tnode = (*tnode) -> left;
		     		      
		  (*tnode) -> left -> level = (*tnode) -> level + 1;
		      
		  for(i = 0; i < 4; i++)
		    {
		      (*tnode) -> left -> bounding_box[i] = bounding_box[i];
		      (*tnode) -> left -> lower_bound[i] = (*tnode) -> lower_bound[i]; /*Lower bound in hypercube (temporary value)*/
		      (*tnode) -> left -> upper_bound[i] = (*tnode) -> upper_bound[i]; /*Upper bound in hypercube*/
		    }

		  /*The upper bound has changed - now the midpoint of this subregion on this axis*/
		      
		  (*tnode) -> left -> upper_bound[((*tnode) -> level)%4] = 
		    0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%4]) + ((*tnode) -> upper_bound[((*tnode) -> level)%4]));     
		      		      
		  (*tnode) -> left -> left = NULL; (*tnode) -> left -> right = NULL;

		  (*tnode) -> left -> parent = (*tnode); (*tnode) -> left -> poly_num = poly_num;

		  return(NOERROR);
		}
	      else return(ERROR);
	    }
	  else /*Node already exists*/
	    {
	      if(add_to_face_ADT(Someface, bounding_box, root_2D_bbox, &((*tnode) -> left), poly_num) == ERROR)
		return(ERROR);
	      else return(NOERROR);
	    }	  
	}
      else 
	{ 
	  if((*tnode) -> right == NULL) 
	    {
	      (*tnode) -> right = malloc(sizeof(struct body2D_tnode));

	      if((*tnode) -> right != NULL)
		{
		  (*tnode) -> right -> Poly = Someface;
		      
		  (*tnode) -> right -> Poly -> tnode = (*tnode) -> right;
		     
		  (*tnode) -> right -> level = (*tnode) -> level + 1;
		      		      	      
		  for(i = 0; i < 4; i++)
		    {
		      (*tnode) -> right -> bounding_box[i] = bounding_box[i];
		      (*tnode) -> right -> lower_bound[i] = (*tnode) -> lower_bound[i]; /*Lower bound in hypercube*/
		      (*tnode) -> right -> upper_bound[i] = (*tnode) -> upper_bound[i]; /*Upper bound in hypercube*/
		    }
		      
		  /*The lower bound has changed - now the midpoint of this subregion on this axis*/

		  (*tnode) -> right -> lower_bound[((*tnode) -> level)%4] = 
		    0.5*(((*tnode) -> lower_bound[((*tnode) -> level)%4]) + ((*tnode) -> upper_bound[((*tnode) -> level)%4]));

		  (*tnode) -> right -> left = NULL; (*tnode) -> right -> right = NULL;

		  (*tnode) -> right -> parent = (*tnode); (*tnode) -> right -> poly_num = poly_num;

		  return(NOERROR);		    
		}
	      else return(ERROR);
	    }
	  else /*Node already exists*/
	    {
	      if(add_to_face_ADT(Someface, bounding_box, root_2D_bbox, &((*tnode) -> right), poly_num) == ERROR)	      
		return(ERROR);
	      else return(NOERROR);
	    }	  
	}
    }
}

/*------------------------------------------------------------------*/
