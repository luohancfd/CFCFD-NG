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

/**\file Data structures and function prototypes for geometry interrogation functions*/

typedef struct body_tnode * Body_tnode;
typedef struct body2D_tnode * Body2D_tnode;
typedef struct edge1D_tnode * Edge1D_tnode;
typedef struct body * Body;
typedef struct body2D * Body2D;
typedef struct edge3D * Edge3D;
typedef struct list_bbox * List_bbox; 

#define PI 3.1415926535897932384626433832795028841972
#define TOL 1e-13 /*Main tolerance of geometry engine (generally related to floating point computations i.e. acos)*/
#define TOL2 0.2 /*Tolerance related to computation of point location relative to polygonal plane (can be more generous here)*/

#define NOTINCLUDE 0

#define TRUE 1
#define FALSE 0
#define NOERROR 0
#define ERROR -1
#define OUT_POLY 2 
#define IN_POLY 3
#define ON_POLY_EDGE 4 /**< Point is on edge of polygon*/
#define ON_POLY_VTX 5
#define OUT_POLYH 6
#define ON_POLYH 7
#define IN_POLYH 8

#define ABS(A) (((A) < 0.0) ? (-(A)): (A)) 
#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) > (B)) ? (A) : (B))
#define SGN(A) (((A) == 0) ? 0 : ((((double) (A))/(ABS((double) (A))))))
#define SQR(A) ((A)*(A))
#define SQRT(A) (pow(((double) (A)), 0.5))
#define LENGTH2(U) (SQRT(SQR(U[0]) + SQR(U[1])))

/*Data structures*/
/*------------------------------------------------------------------*/

struct body_tnode /**< \b ADT of bodies (and their bounding boxes)*/
{
  Body Solid; /**< Store solid body*/
  int level; /**< Level of node on ADT*/
    
  double bounding_box[6]; /**< Bounding box in hyperspace*/
  double plane_bbox[3][4]; /**< 4D Projected bounding box on plane (row-wise are the 3 directions)*/
  
  double lower_bound[6];
  double upper_bound[6];
  /**< Partitioned region corresponding to body in hypercube*/

  struct body_tnode *left; struct body_tnode *right;
  
  int body_num; /*body_num th body to be put on ADT*/
  struct body_tnode *parent;
};

struct body2D_tnode /**< ADT of 2D bodies (and their bounding boxes)*/
{
  Body2D Poly;
  int level;

  double bounding_box[4];
  double lower_bound[4];
  double upper_bound[4];

  struct body2D_tnode *left; struct body2D_tnode *right;

  int poly_num;
  struct body2D_tnode *parent;
};

struct edge1D_tnode /**< ADT of projected polygon edges (on x axis)*/
{
  Edge3D Edge;
  int level;

  double bounding_box[2];
  double lower_bound[2];
  double upper_bound[2];

  struct edge1D_tnode *left; struct edge1D_tnode *right;

  int edge_num;
  struct edge1D_tnode *parent;
};

struct body /**< \b Store geometric info on body*/
{
  Body_tnode tnode; /**< Pointer to corresponding node on body ADT*/
  Body2D_tnode Faces; /**< Root node of ADT of projected body facets*/
  int num_points;
  int num_faces;

  double xy_bound[4]; /**< Bounds on XY plane that projected face polygons occupy (in 4D space)*/

  double **points; /**< Array of all points (verticies) defining polyhedron; column-wise are co-ordinates*/
  int **faces; /**< Matrix of faces (1 per row) referencing points[0] array; faces[i][0] start point, faces[i][n] end point;
		  faces defined in right-hand system with outward pointing normal (though normal not needed)*/
  int *face_size; /**< Each element of face_size array denotes how many verticies for ith face in faces array*/
};

struct body2D /**< 2D body (can exist in 3D space but has no 'thickness')*/
{
  Body2D_tnode tnode;
  Edge1D_tnode Edges; /**< Root node of ADT of projected edges*/
  int num_points; /**< Non-degenerate polygons have as many points as edges*/
  
  double x_bound[2]; /**< Bounds over x-axis which this body occupies (in 2D space)*/

  double **points; /**< Array of all points in polygon (row-wise); column-wise are co-ordinates*/
  int **edges; /**< Matrix of edges (1 per row) referencing points[] array; edges[i][0] start point, edges[i][1] end point;
		  edges defined clockwise/counter-clockwise*/
  double norm[3]; /**< Components of surface normal to polygon*/

#if NOTINCLUDE
  int poly_num;
#endif
};

struct edge3D /**< Edge of polygon*/
{
  Edge1D_tnode tnode;
  Body2D Poly; /**< Belongs on this polygon*/
  double start[3];
  double end[3];
};

struct list_bbox /**< \b List of body bounding boxes (intersecting given region)*/
{
  Body_tnode tnode; /**< Stores address of bounding box*/
  struct list_bbox *next;
};

/*Geometry engine*/
/*------------------------------------------------------------------*/

short int polyhio(double[], Body);

short int polyio(double[], Body2D, double *);

short int get_link_num(Edge3D, double[], int *);

double compute_side_loc(double[], Body2D);

short int precompute_normal(Body2D);

double get_vtx_angle(Body2D, int, int);

short int collinear_test(double[], double[], double[]);

/*------------------------------------------------------------------*/
