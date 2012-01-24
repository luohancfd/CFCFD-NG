#include <stdio.h>
#include <stdlib.h>
#include "gio_kernel.h"
#include "gio_adts.h"
#include "gio_lists.h"

#define RECTANGLE 0
#define LINE 0
#define DIAMOND 1
#define HEXAGON 0
#define PENTAGON 0

void get_point(double []);
short int make_polygon(Body2D);
void traverse_edges(Edge1D_tnode);

main()
{
  Body2D Poly;
  Edge3D Edge;
  double bbox[2];
  int i;
  double root_bounds[2][1];
  double point[2];

  Poly = malloc(sizeof(struct body2D));
  if(Poly == NULL)
    return(ERROR);

  Poly -> Edges = NULL; /*Must safely initialize*/
  
  /*Define the polygon*/

  get_point(point);
  make_polygon(Poly);
  
  /*Get minimum and maximum bounds (on x-axis) for polygon*/

  root_bounds[0][0] = Poly -> points[0][0];
  root_bounds[1][0] = root_bounds[0];

  for(i = 0; i < Poly -> num_points; i++)
    {
      if(Poly -> points[i][0] < root_bounds[0][0])
	root_bounds[0][0] = Poly -> points[i][0];
      else if(Poly -> points[i][0] > root_bounds[1][0])
	root_bounds[1][0] = Poly -> points[i][0];
    }

  Poly -> x_bound[0] = root_bounds[0][0]; 
  Poly -> x_bound[1] = root_bounds[1][0];

  /*Add edges to ADT of edges*/

  for(i = 0; i < Poly -> num_edges; i++)
    {
      Edge = malloc(sizeof(struct edge3D));
      if(Edge == NULL)
	return(ERROR);

      Edge -> start[0] = Poly -> points[Poly -> edges[i][0]][0]; /*Edge matrix only references points which constitute its start/end*/
      Edge -> start[1] = Poly -> points[Poly -> edges[i][0]][1];
      Edge -> end[0] = Poly -> points[Poly -> edges[i][1]][0];
      Edge -> end[1] = Poly -> points[Poly -> edges[i][1]][1];

      bbox[0] = MIN(Edge -> start[0],Edge -> end[0]); /*Edges might be defined clockwise/counter-clockwise*/
      bbox[1] = MAX(Edge -> start[0],Edge -> end[0]);
      
      if(add_to_edge_ADT(Edge, bbox, root_bounds, &(Poly -> Edges), i+1) == ERROR)
	{
	  printf("Add_to_edge_ADT(): Error.  Exiting ...\n");
	  return(ERROR);
	}

      Edge = NULL;
    }

  /*Now make query*/

  i = polyio(point, Poly);
  if(i == IN_POLY)
    printf("INSIDE POLYGON\n");
  else if(i == OUT_POLY)
    printf("OUTSIDE POLYGON\n");
}

void get_point(double point[]) /*Define points here for the moment*/
{
  point[0] = -0.5  ;
  point[1] = 0.4999  ;

  printf("The point I'm queryin' is (%g,%g)\n",point[0],point[1]);
}

short int make_polygon(Body2D Poly) /*Define polygon here for the moment*/
{
  int i;

#if RECTANGLE
  Poly -> num_points = 4;
  Poly -> num_edges = 4;
#elif LINE
  Poly -> num_points = 2; /*Try if degenerate line is OK*/
  Poly -> num_edges = 4;
#elif DIAMOND
  Poly -> num_points = 4;
  Poly -> num_edges = 4;
#elif HEXAGON
  Poly -> num_points = 6;
  Poly -> num_edges = 6;
#elif PENTAGON
  Poly -> num_points = 5;
  Poly -> num_edges = 5;
#endif

  Poly -> points = malloc(sizeof(double *)*(Poly -> num_points));
  if(Poly -> points == NULL)
    return(ERROR);

  for(i = 0; i < Poly -> num_points; i++) 
    {
      Poly -> points[i] = malloc(sizeof(double)*3); /*Point in 3D space, but will only project onto XY plane anyway*/
      if(Poly -> points[i] == NULL)
	return(ERROR);
    }

  Poly -> edges = malloc(sizeof(int *)*(Poly -> num_edges));
  if(Poly -> edges == NULL)
    return(ERROR);

  for(i = 0; i < Poly -> num_edges; i++)
    {
      Poly -> edges[i] = malloc(sizeof(int)*2);
      if(Poly -> edges[i] == NULL)
	return(ERROR);
    }

  /*Define points and edges here*/

#if RECTANGLE
  Poly -> points[0][0] = 0; Poly -> points[0][1] = 0; /*1st point, etc.*/
  Poly -> points[1][0] = 2; Poly -> points[1][1] = 0;
  Poly -> points[2][0] = 2; Poly -> points[2][1] = 1;
  Poly -> points[3][0] = 0; Poly -> points[3][1] = 1;

  Poly -> edges[0][0] = 0; Poly -> edges[0][1] = 1; /*1st edge (start and end points), etc ...*/
  Poly -> edges[1][0] = 1; Poly -> edges[1][1] = 2;
  Poly -> edges[2][0] = 2; Poly -> edges[2][1] = 3;
  Poly -> edges[3][0] = 3; Poly -> edges[3][1] = 0;

#elif LINE
  Poly -> points[0][0] = 0; Poly -> points[0][1] = 0;
  Poly -> points[1][0] = 1; Poly -> points[1][1] = 1;
  /*Poly -> points[2][0] = 2; Poly -> points[2][1] = 2;*/

  Poly -> edges[0][0] = 0; Poly -> edges[0][1] = 1;
  Poly -> edges[1][0] = 1; Poly -> edges[1][1] = 1;
  Poly -> edges[2][0] = 1; Poly -> edges[2][1] = 0;
  Poly -> edges[3][0] = 0; Poly -> edges[3][1] = 0;

#elif DIAMOND /*A clockwise diamond*/
  Poly -> points[0][0] = 0; Poly -> points[0][1] = 0; Poly -> points[0][2] = -0.234; /*An irrelevant height*/
  Poly -> points[1][0] = 1; Poly -> points[1][1] = 1; Poly -> points[1][2] = 1.234;
  Poly -> points[2][0] = 0; Poly -> points[2][1] = 2; Poly -> points[2][2] = -1234;
  Poly -> points[3][0] = -1; Poly -> points[3][1] = 1; Poly -> points[3][2] = 5;

  Poly -> edges[0][0] = 0; Poly -> edges[0][1] = 3;
  Poly -> edges[1][0] = 3; Poly -> edges[1][1] = 2;
  Poly -> edges[2][0] = 2; Poly -> edges[2][1] = 1;
  Poly -> edges[3][0] = 1; Poly -> edges[3][1] = 0;

#elif HEXAGON
  Poly -> points[0][0] = 0; Poly -> points[0][1] = 0;
  Poly -> points[1][0] = 1; Poly -> points[1][1] = 1;
  Poly -> points[2][0] = 2; Poly -> points[2][1] = 1;
  Poly -> points[3][0] = 3; Poly -> points[3][1] = 0;
  Poly -> points[4][0] = 2; Poly -> points[4][1] = -1;
  Poly -> points[5][0] = 1; Poly -> points[5][1] = -1;

  Poly -> edges[0][0] = 0; Poly -> edges[0][1] = 5;
  Poly -> edges[1][0] = 5; Poly -> edges[1][1] = 4;
  Poly -> edges[2][0] = 4; Poly -> edges[2][1] = 3;
  Poly -> edges[3][0] = 3; Poly -> edges[3][1] = 2;
  Poly -> edges[4][0] = 2; Poly -> edges[4][1] = 1;
  Poly -> edges[5][0] = 1; Poly -> edges[5][1] = 0;

#elif PENTAGON /*Non-convex test*/
  Poly -> points[0][0] = 0; Poly -> points[0][1] = 0;
  Poly -> points[1][0] = -1; Poly -> points[1][1] = -1;
  Poly -> points[2][0] = 1; Poly -> points[2][1] = -1;
  Poly -> points[3][0] = 1; Poly -> points[3][1] = 1;
  Poly -> points[4][0] = -1; Poly -> points[4][1] = 1;
  
  Poly -> edges[0][0] = 0; Poly -> edges[0][1] = 1;
  Poly -> edges[1][0] = 1; Poly -> edges[1][1] = 2;
  Poly -> edges[2][0] = 2; Poly -> edges[2][1] = 3;
  Poly -> edges[3][0] = 3; Poly -> edges[3][1] = 4;
  Poly -> edges[4][0] = 4; Poly -> edges[4][1] = 0;
#endif

  printf("So the points of the polygon in order are ...\n");
  for(i = 0; i < Poly -> num_points; i++)
    printf("(%g,%g)\n",Poly -> points[i][0], Poly -> points[i][1]);

  return(NOERROR);
}

void traverse_edges(Edge1D_tnode tnode)
{
  if(tnode -> parent != NULL) {
    printf("Edge %d, Parent %d, level %d, bound [%g,%g], ",tnode->edge_num,tnode->parent->edge_num,tnode->level,tnode->bounding_box[0],
	   tnode->bounding_box[1]);
    if(tnode == tnode -> parent -> left)
      printf("left, ");
    else if(tnode == tnode -> parent -> right)
      printf("right, ");
    printf("Lwr[%g,%g], Upr[%g,%g]\n",tnode->lower_bound[0],tnode->lower_bound[1],tnode->upper_bound[0],tnode->upper_bound[1]);
  }
  else {
    printf("Edge %d, level %d, bound [%g,%g], ",tnode->edge_num,tnode->level,tnode->bounding_box[0],tnode->bounding_box[1]);
    printf("Lwr[%g,%g], Upr[%g,%g]\n",tnode->lower_bound[0],tnode->lower_bound[1],tnode->upper_bound[0],tnode->upper_bound[1]);
  }

  if(tnode -> left != NULL)
    traverse_edges(tnode->left);
  if(tnode -> right != NULL)
    traverse_edges(tnode->right);
}
