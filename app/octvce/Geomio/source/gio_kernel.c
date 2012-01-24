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


/**\file Source file for all global variables and interrogation engine - hopefully polygons and polyhedrons.  These set of routines
   were written primarily to be used in conjunction with Octree Virtual Cell Embedding (OctVCE).  Point inclusion queries here
   implement the algorithms set forth in 

   1. Linhart, J, "A quick point-in-polyhedron test", Computers and Graphics, vol. 14(3/4), pp. 445-447, 1990

   2. Wu, H, Gong, J, Li, D & Shi, W, "An algebraic algorithm for point inclusion query", Computers and Graphics, 
      vol. 24(4), pp. 517-522, 2000.

   I've also modified the algorithm in [1] slightly to make it more similar to [2]; the geometry engine may thus classify some points
   as being inside the polyhedron though they're on the surface, but this shouldn't matter too much (as far as OctVCE is concerned).

   As long as we don't have highly skewed shapes, the algorithm is robust, in that no point outside the surface will be calculated
   as on/inside it (given machine and numerical tolerance TOL).
   
   The point inclusion query procedures coded here have been significantly preprocessed, sacrificing memory and storage in favour
   of speed, and are thus suitable for fast(er) queries involving multiple test points (as in OctVCE) or polyhedra with many edges.      
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gio_kernel.h"
#include "gio_adts.h"

#define NOTINCLUDE 0

Body_tnode ADTbody = NULL; /**< ADT of bodies (initialized to NULL - if NULL, no bodies)*/
double root_3D_bbox[2][3] = {{-1e20,-1e20,-1e20}, {1e20,1e20,1e20}}; /*Bounding box of root domain (of body ADT)*/
int body_count = 0; /*No. bodies read so far*/

/*------------------------------------------------------------------*/

/**\brief Test if 3D point inside arbitrary polyhedron*/

short int polyhio(double point[], Body Solid)
{
  double s = 0; /*Weight for particular polygon (s == 0, point is outside)*/
  double lower_bound[4];
  double upper_bound[4];
  short int i;

  /*Get bounds on degenerate 4D hyperrectangle which the point occupies*/

  lower_bound[0] = Solid -> xy_bound[0];
  lower_bound[1] = Solid -> xy_bound[1];
  lower_bound[2] = point[0];
  lower_bound[3] = point[1];

  upper_bound[0] = point[0];
  upper_bound[1] = point[1];
  upper_bound[2] = Solid -> xy_bound[2];
  upper_bound[3] = Solid -> xy_bound[3]; 

  i = traverse_face_ADT(point, lower_bound, upper_bound, Solid -> Faces, &s); 
  
  if(i == ON_POLYH)
    return(i);
  else if((s >= -TOL) && (s <= TOL)) /*With floating point computations must give some tolerance*/
    return(OUT_POLYH);
  else if(s >= -1) /*If point lies on vertical edge, it'll pass through 2 verticies with varying inner angles - sum won't be < -1*/
    return(IN_POLYH);
  else return(ERROR);
}

/*------------------------------------------------------------------*/

/**\brief Test if 2D point inside arbitrary polygon*/

short int polyio(double point[], Body2D Poly, double *vtx_angle)
{
  int link_num = 0; /*Link number of point to polygon (projected onto XY plane)*/
  double lower_bound[2]; /*Point represented as degenerate area in 2 dimensional space*/
  double upper_bound[2];
  short int i;

  lower_bound[0] = Poly -> x_bound[0]; 
  lower_bound[1] = point[0];
  upper_bound[0] = point[0]; 
  upper_bound[1] = Poly -> x_bound[1]; /*Only interested in ranges over x-coordinate*/

  i = traverse_edge_ADT(point, lower_bound, upper_bound, Poly -> Edges, &link_num, vtx_angle); 
  
  if((i == ON_POLY_VTX) || (i == ON_POLY_EDGE))
    return(i);
  else if(link_num == 0)
    return(OUT_POLY);
  else if(ABS(link_num) == 4)
    return(IN_POLY);
  else if(ABS(link_num) == 2)
    return(ON_POLY_EDGE); /*If it were on the vertex then i == ON_POLY_VTX already*/
  else return(ERROR); /*A degenerate polygon*/
}

/*------------------------------------------------------------------*/

/**\brief Calculate linking number of a point to an edge segment (in 2D, on XY plane)*/

short int get_link_num(Edge3D Edge, double point[], int *link_num)
{
  int cross = 2;
  double nDist;

  if((Edge->start[0] != Edge->end[0]) && 
     (((Edge->start[0] <= point[0]) && (Edge->end[0] >= point[0])) || ((Edge->start[0] >= point[0]) && (Edge->end[0] <= point[0]))))
    { /*A crossing exists and isn't degenerate (non-crossing)*/
      
      if((Edge -> start[0] == point[0]) || (Edge -> end[0] == point[0]))
	cross = 1; /*Semicrossing*/
      nDist = (point[1] - Edge->start[1])*(Edge->end[0] - Edge->start[0]) - (point[0] - Edge->start[0])*(Edge->end[1] - Edge->start[1]); 

      /*Shortest distance from point to line (in 2D)*/

      if(nDist == 0) 
	{
	  if(cross == 1)
	    return(ON_POLY_VTX);
	  else return(ON_POLY_EDGE);
	}
      else if(nDist > 0)
	*link_num += cross; /*Positive crossing*/
      else *link_num -= cross; /*Negative crossing*/      
    }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Determine if a point is on upside/downside of a plane defined by a polygon*/

double compute_side_loc(double point[], Body2D Poly)
{
  short int i;
  double p[3];
  double sign;

  for(i = 0; i < 3; i++)
    p[i] = point[i] - (Poly -> points[Poly->edges[0][0]][i]); /*For simplicity chose 0th point of polygon as point on plane*/

  sign = SGN(p[0]*(Poly->norm[0]) + p[1]*(Poly->norm[1]) + p[2]*(Poly->norm[2]));

  if((sign >= -TOL2) && (sign <= TOL2))
    sign = 0; 
  else if((sign >= (-1.0-TOL2)) && (sign <= (-1.0+TOL2)))
    sign = -1.0; 
  else if ((sign >= (1.0-TOL2)) && (sign <= (1.0+TOL2))) 
    sign = 1.0;

  return(sign); /*Sign of scalar triple product computed*/
}

/*------------------------------------------------------------------*/

/**\brief Precompute normal to surface of polygonal facet for later scalar triple product calculations*/

short int precompute_normal(Body2D Poly)
{
  double u[3], v[3];
  int i;
  short int j = TRUE;

  for(i = 1; i <= (Poly -> num_points - 2); i++) /*Take 1st two points and find a third non-collinear one*/
    {
      j = collinear_test(Poly -> points[Poly->edges[0][0]], Poly -> points[Poly->edges[0][1]], Poly -> points[Poly->edges[i][1]]);
      if(j == FALSE)
	break;
    }

  if(j == TRUE) /*Couldn't find any non-collinear point - this is a degenerate polygon*/
    return(ERROR);

  for(j = 0; j < 3; j++) /*Form vectors to cross*/
    {
      u[j] = (Poly -> points[Poly->edges[0][1]][j]) - (Poly -> points[Poly->edges[0][0]][j]);
      v[j] = (Poly -> points[Poly->edges[i][1]][j]) - (Poly -> points[Poly->edges[i][0]][j]); 
    }

  Poly -> norm[0] = u[1]*v[2] - v[1]*u[2];
  Poly -> norm[1] = v[0]*u[2] - u[0]*v[2];
  Poly -> norm[2] = u[0]*v[1] - v[0]*u[1]; /*Normal computed*/

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Get inner vertex angle of projected polygon (on xy plane) between edges*/

double get_vtx_angle(Body2D Poly, int e0, int e1)
{
  short int i;
  double angle = 0;
  double norm, u[2], v[2];

  for(i = 0; i < 2; i++)
    {
      u[i] = (Poly -> points[Poly -> edges[e0][1]][i]) - (Poly -> points[Poly -> edges[e0][0]][i]); /*Vector of 1st edge*/
      v[i] = (Poly -> points[Poly -> edges[e1][1]][i]) - (Poly -> points[Poly -> edges[e1][0]][i]); /*Vector of 2nd edge*/
    }  

  norm = u[0]*v[1] - u[1]*v[0]; /*Z-component of normal vector to u and v*/
  
  if(SGN(Poly -> norm[2]) == 0) 
    return(0); /*Polygon plane orthogonal to xy plane - angle is 0 (traverse_face_ADT would see the point lies on the polygon)*/
  
  if(SGN(Poly -> norm[2]) < 0) /*Polygon is oriented 'downward', so edges on xy plane are traversed clockwise*/
    angle = ABS(PI + SGN(norm)*acos((u[0]*v[0]+u[1]*v[1])/((LENGTH2(u))*(LENGTH2(v))))); /*Use dot product as sine is ambiguous*/
  else angle = ABS(-PI + SGN(norm)*acos((u[0]*v[0]+u[1]*v[1])/((LENGTH2(u))*(LENGTH2(v)))));

  return(angle);
}

/*------------------------------------------------------------------*/

/**\brief Test for collinearity of 3 consecutive points in 3D*/

short int collinear_test(double p0[], double p1[], double p2[])
{
  /*It is assumed that no points are identical.  Test collinearity based
    on dot product (magnitude of dot product would equal vector lengths multiplied by each other)*/

  double u[3], v[3];

  u[0] = p1[0] - p0[0];
  u[1] = p1[1] - p0[1];
  u[2] = p1[2] - p0[2];

  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  if(SQR(u[0]*v[0] + u[1]*v[1] + u[2]*v[2]) == (SQR(u[0])+SQR(u[1])+SQR(u[2]))*(SQR(v[0])+SQR(v[1])+SQR(v[2])))
    return(TRUE);
  
  return(FALSE); /*Non-collinear*/
}

/*------------------------------------------------------------------*/


