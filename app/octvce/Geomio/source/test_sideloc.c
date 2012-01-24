#include <stdio.h>
#include <stdlib.h>
#include "gio_kernel.h"

#define ZPLANE 0
#define SLANTED_PLANE 1

void get_point(double []);
short int make_polygon(Body2D);

main()
{
  short int i;
  double point[3];
  Body2D Poly = malloc(sizeof(struct body2D));

  if(Poly == NULL)
    return(ERROR);

  /*Now define the point and polygon*/
  make_polygon(Poly);
  get_point(point);
  
  if(precompute_normal(Poly) == ERROR)
    return(ERROR);

  printf("Dets are %g %g %g\n",Poly->norm[0],Poly->norm[1],Poly->norm[2]);

  i = compute_side_loc(point, Poly);

  if(i == 0)
    printf("ON POLYGON\n");
  else if(i == -1)
    printf("BENEATH POLYGON\n");
  else if(i == 1)
    printf("ABOVE POLYGON\n");
}

void get_point(double point[]) /*Define points here for the moment*/
{
  point[0] = -10  ;
  point[1] = 20  ;
  point[2] = 19.999999   ;

  printf("The point I'm queryin' is (%g,%g, %g)\n",point[0],point[1],point[2]);
}

short int make_polygon(Body2D Poly) /*Define polygon here for the moment*/
{
  int i;

#if ZPLANE
  Poly -> num_points = 5;
#elif SLANTED_PLANE
  Poly -> num_points = 6;
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

  Poly -> edges = malloc(sizeof(int *)*(Poly -> num_points));
  if(Poly -> edges == NULL)
    return(ERROR);

  for(i = 0; i < Poly -> num_points; i++)
    {
      Poly -> edges[i] = malloc(sizeof(int)*2);
      if(Poly -> edges[i] == NULL)
	return(ERROR);
    }

  /*Define points and edges here*/

#if ZPLANE
  /*Rectangle with 2 collinear edges (on z-plane)*/
  Poly -> points[0][0] = 0; Poly -> points[0][1] = 0; Poly -> points[0][2] = 0;
  Poly -> points[1][0] = 1; Poly -> points[1][1] = 0; Poly -> points[1][2] = 0;
  Poly -> points[2][0] = 2; Poly -> points[2][1] = 0; Poly -> points[2][2] = 0;
  Poly -> points[3][0] = 2; Poly -> points[3][1] = 1; Poly -> points[3][2] = 0;
  Poly -> points[4][0] = 0; Poly -> points[4][1] = 1; Poly -> points[4][2] = 0;

  Poly -> edges[0][0] = 0; Poly -> edges[0][1] = 1;
  Poly -> edges[1][0] = 1; Poly -> edges[1][1] = 2;
  Poly -> edges[2][0] = 2; Poly -> edges[2][1] = 3;
  Poly -> edges[3][0] = 3; Poly -> edges[3][1] = 4;
  Poly -> edges[4][0] = 4; Poly -> edges[4][1] = 0;

#elif SLANTED_PLANE
  /*Slanted plane (on YZ plane)*/
  Poly -> points[0][0] = 0; Poly -> points[0][1] = 0; Poly -> points[0][2] = 0;
  Poly -> points[1][0] = 0.5; Poly -> points[1][1] = 0; Poly -> points[1][2] = 0;
  Poly -> points[2][0] = 1; Poly -> points[2][1] = 0; Poly -> points[2][2] = 0;
  Poly -> points[3][0] = 2; Poly -> points[3][1] = 0; Poly -> points[3][2] = 0;
  Poly -> points[4][0] = 2; Poly -> points[4][1] = 1; Poly -> points[4][2] = 1;
  Poly -> points[5][0] = 0; Poly -> points[5][1] = 1; Poly -> points[5][2] = 1;

  Poly -> edges[0][0] = 0; Poly -> edges[0][1] = 1;
  Poly -> edges[1][0] = 1; Poly -> edges[1][1] = 2;
  Poly -> edges[2][0] = 2; Poly -> edges[2][1] = 3;
  Poly -> edges[3][0] = 3; Poly -> edges[3][1] = 4;
  Poly -> edges[4][0] = 4; Poly -> edges[4][1] = 5;
  Poly -> edges[5][0] = 5; Poly -> edges[5][1] = 0;
#endif

  printf("So the points of the polygon in order are ...\n");
  for(i = 0; i < Poly -> num_points; i++)
    printf("(%g,%g,%g)\n",Poly -> points[i][0], Poly -> points[i][1], Poly -> points[i][2]);

  return(NOERROR);
}
