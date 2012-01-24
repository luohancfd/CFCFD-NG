#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gio_kernel.h"
#include "gio_adts.h"
#include "gio_io.h"

#define RECTANGULAR_PRISM 0
#define TRIANGULAR_PRISM 0
#define DEGENERATE_TETRAHEDRON 0
#define CONCAVE_SHAPE 0
#define DEGENERATE_RECT_PRISM 1

void get_point(double []);
void traverse_faces(Body2D_tnode);
short int make_polyh(Body);

extern Body_tnode ADTbody;

main()
{
  Body Solid;
  double point[3];
  int i, j;

  if(read_body_files("bodies.pvtp.xml") == ERROR) {
    printf("Unable to successfully read the body files\n");
    return(ERROR);    
  }

  /*So we've already built ADTbody (for the moment just 1 body).  Now can interrogate it*/

  Solid = ADTbody -> Solid;
  printf("The solid has %d points and %d faces\n",Solid->num_points,Solid->num_faces);
  printf("Points: \n");
  for(i = 0; i < Solid->num_points; i++)
    printf("%g %g %g\n",Solid->points[i][0], Solid->points[i][1], Solid->points[i][2]);
  printf("Faces: \n");
  for(i = 0; i < Solid->num_faces; i++)
    {
      for(j = 0; j < Solid->face_size[i]; j++)
	printf("%d ", Solid->faces[i][j]);
      printf("\n");
    }
  
  get_point(point);
  
#if 0
  Solid = malloc(sizeof(struct body));
  if(Solid == NULL)
    return(ERROR);

  /*Define polyhedron*/
  get_point(point);
  make_polyh(Solid);
  
  /*Build ADT*/
  if(build_face_ADT(Solid, TRUE) == ERROR)
    return(ERROR);
#endif
  
  i = polyhio(point, Solid); 
  if(i == IN_POLYH)
    printf("INSIDE\n");
  else if(i == OUT_POLYH)
    printf("OUTSIDE\n");
  else if(i == ON_POLYH)
    printf("ON\n");
  else printf("ERROR\n");
  
  return(NOERROR);
}

void get_point(double point[]) /*Define points here for the moment*/
{
  point[0] = -0.2 ;
  point[1] = -0.15 ;
  point[2] =  0.21000001 ;

  printf("The point I'm queryin' is (%g,%g,%g)\n",point[0],point[1],point[2]);
}

short int make_polyh(Body Solid)
{
  int i, j;

#if RECTANGULAR_PRISM
  Solid -> num_points = 8;
  Solid -> num_faces = 6;

#elif TRIANGULAR_PRISM
  Solid -> num_points = 6;
  Solid -> num_faces = 5;

#elif DEGENERATE_TETRAHEDRON
  Solid -> num_points = 5;
  Solid -> num_faces = 5;

#elif CONCAVE_SHAPE /*8-faced shaped body looking like a pyramid inside a pyramid*/
  Solid -> num_points = 6;
  Solid -> num_faces = 8;

#elif DEGENERATE_RECT_PRISM
  Solid -> num_points = 13;
  Solid -> num_faces = 9;
#endif

  Solid -> points = malloc(sizeof(double *)*(Solid -> num_points));
  if(Solid -> points == NULL)
    return(ERROR);

  for(i = 0; i < Solid -> num_points; i++) 
    {
      Solid -> points[i] = malloc(sizeof(double)*3);
      if(Solid -> points[i] == NULL)
	return(ERROR);
    }

  Solid -> faces = malloc(sizeof(int *)*(Solid -> num_faces));
  if(Solid -> faces == NULL)
    return(ERROR);

  Solid -> face_size = malloc(sizeof(int)*(Solid -> num_faces));
  if(Solid -> face_size == NULL)
    return(ERROR);

  /*Define points and faces here - must allocate size of face*/

#if RECTANGULAR_PRISM
  Solid -> points[0][0] = 0; Solid -> points[0][1] = 0; Solid -> points[0][2] = 0;
  Solid -> points[1][0] = 0; Solid -> points[1][1] = 1; Solid -> points[1][2] = 0;
  Solid -> points[2][0] = 1; Solid -> points[2][1] = 1; Solid -> points[2][2] = 0;
  Solid -> points[3][0] = 1; Solid -> points[3][1] = 0; Solid -> points[3][2] = 0;
  Solid -> points[4][0] = 0; Solid -> points[4][1] = 0; Solid -> points[4][2] = 2;
  Solid -> points[5][0] = 0; Solid -> points[5][1] = 1; Solid -> points[5][2] = 2;
  Solid -> points[6][0] = 1; Solid -> points[6][1] = 1; Solid -> points[6][2] = 2;
  Solid -> points[7][0] = 1; Solid -> points[7][1] = 0; Solid -> points[7][2] = 2;

  for(i = 0; i < Solid -> num_faces; i++)
    Solid -> face_size[i] = 4;

  for(i = 0; i < Solid -> num_faces; i++) {
    Solid -> faces[i] = malloc(sizeof(int)*(Solid -> face_size[i]));
    if(Solid -> faces[i] == NULL)
      return(ERROR);
  }

  Solid -> faces[0][0] = 6; Solid -> faces[0][1] = 5; Solid -> faces[0][2] = 4; Solid -> faces[0][3] = 7;
  Solid -> faces[1][0] = 2; Solid -> faces[1][1] = 6; Solid -> faces[1][2] = 7; Solid -> faces[1][3] = 3;
  Solid -> faces[2][0] = 0; Solid -> faces[2][1] = 4; Solid -> faces[2][2] = 5; Solid -> faces[2][3] = 1;
  Solid -> faces[3][0] = 1; Solid -> faces[3][1] = 5; Solid -> faces[3][2] = 6; Solid -> faces[3][3] = 2;
  Solid -> faces[4][0] = 0; Solid -> faces[4][1] = 3; Solid -> faces[4][2] = 7; Solid -> faces[4][3] = 4;
  Solid -> faces[5][0] = 1; Solid -> faces[5][1] = 2; Solid -> faces[5][2] = 3; Solid -> faces[5][3] = 0;  

#elif TRIANGULAR_PRISM
  Solid -> points[0][0] = 0; Solid -> points[0][1] = 0; Solid -> points[0][2] = 0;
  Solid -> points[1][0] = 0; Solid -> points[1][1] = 1; Solid -> points[1][2] = 0;
  Solid -> points[2][0] = 1; Solid -> points[2][1] = 0; Solid -> points[2][2] = 0;
  Solid -> points[3][0] = 0; Solid -> points[3][1] = 0; Solid -> points[3][2] = 1;
  Solid -> points[4][0] = 0; Solid -> points[4][1] = 1; Solid -> points[4][2] = 1;
  Solid -> points[5][0] = 1; Solid -> points[5][1] = 0; Solid -> points[5][2] = 1;

  Solid -> face_size[0] = 3;
  Solid -> face_size[1] = 4;
  Solid -> face_size[2] = 4;
  Solid -> face_size[3] = 4;
  Solid -> face_size[4] = 3;

  for(i = 0; i < Solid -> num_faces; i++) {
    /*Allocate each face row*/
    Solid -> faces[i] = malloc(sizeof(int)*(Solid -> face_size[i]));
    if(Solid -> faces[i] == NULL)
      return(ERROR);
  }

  Solid -> faces[0][0] = 0; Solid -> faces[0][1] = 1; Solid -> faces[0][2] = 2; 
  Solid -> faces[1][0] = 0; Solid -> faces[1][1] = 2; Solid -> faces[1][2] = 5; Solid -> faces[1][3] = 3;
  Solid -> faces[2][0] = 1; Solid -> faces[2][1] = 4; Solid -> faces[2][2] = 5; Solid -> faces[2][3] = 2;
  Solid -> faces[3][0] = 3; Solid -> faces[3][1] = 4; Solid -> faces[3][2] = 1; Solid -> faces[3][3] = 0;
  Solid -> faces[4][0] = 4; Solid -> faces[4][1] = 3; Solid -> faces[4][2] = 5;

#elif DEGENERATE_TETRAHEDRON
  /*A tetrahedron with 1 more edge than necessary*/

  Solid -> points[0][0] = 1; Solid -> points[0][1] = 0; Solid -> points[0][2] = 0; 
  Solid -> points[1][0] = 2; Solid -> points[1][1] = 0; Solid -> points[1][2] = 0;
  Solid -> points[2][0] = 3; Solid -> points[2][1] = 0; Solid -> points[2][2] = 0;
  Solid -> points[3][0] = 2; Solid -> points[3][1] = SQRT(3); Solid -> points[3][2] = 0;
  Solid -> points[4][0] = 2; Solid -> points[4][1] = SQRT(3)/2; Solid -> points[4][2] = SQRT(3);

  Solid -> face_size[0] = 3;
  Solid -> face_size[1] = 3;
  Solid -> face_size[2] = 3;
  Solid -> face_size[3] = 3;
  Solid -> face_size[4] = 4;
  
  for(i = 0; i < Solid -> num_faces; i++) 
    {
      Solid -> faces[i] = malloc(sizeof(int)*(Solid -> face_size[i]));
      if(Solid -> faces[i] == NULL)
	return(ERROR);
    }

  Solid -> faces[0][0] = 0; Solid -> faces[0][1] = 1; Solid -> faces[0][2] = 4;
  Solid -> faces[1][0] = 1; Solid -> faces[1][1] = 2; Solid -> faces[1][2] = 4;
  Solid -> faces[2][0] = 4; Solid -> faces[2][1] = 2; Solid -> faces[2][2] = 3;
  Solid -> faces[3][0] = 3; Solid -> faces[3][1] = 0; Solid -> faces[3][2] = 4;
  Solid -> faces[4][0] = 3; Solid -> faces[4][1] = 2; Solid -> faces[4][2] = 1; Solid -> faces[4][3] = 0;

#elif CONCAVE_SHAPE

  Solid -> points[0][0] = 0; Solid -> points[0][1] = 0; Solid -> points[0][2] = 0;
  Solid -> points[1][0] = 0; Solid -> points[1][1] = 1; Solid -> points[1][2] = 0;
  Solid -> points[2][0] = 1; Solid -> points[2][1] = 1; Solid -> points[2][2] = 0;
  Solid -> points[3][0] = 1; Solid -> points[3][1] = 0; Solid -> points[3][2] = 0;
  Solid -> points[4][0] = 0.5; Solid -> points[4][1] = 0.5; Solid -> points[4][2] = 1;
  Solid -> points[5][0] = 0.5; Solid -> points[5][1] = 0.5; Solid -> points[5][2] = 3;

  for(i = 0; i < Solid -> num_faces; i++)
    Solid -> face_size[i] = 3;

  for(i = 0; i < Solid -> num_faces; i++) {
    Solid -> faces[i] = malloc(sizeof(int)*(Solid -> face_size[i]));
    if(Solid -> faces[i] == NULL)
      return(ERROR);
  }

  /*Deal with bottom faces first*/
  Solid -> faces[0][0] = 1; Solid -> faces[0][1] = 4; Solid -> faces[0][2] = 0;
  Solid -> faces[1][0] = 0; Solid -> faces[1][1] = 4; Solid -> faces[1][2] = 3;
  Solid -> faces[2][0] = 4; Solid -> faces[2][1] = 2; Solid -> faces[2][2] = 3;
  Solid -> faces[3][0] = 1; Solid -> faces[3][1] = 2; Solid -> faces[3][2] = 4;

  /*Now look at top faces*/
  Solid -> faces[4][0] = 1; Solid -> faces[4][1] = 0; Solid -> faces[4][2] = 5;
  Solid -> faces[5][0] = 0; Solid -> faces[5][1] = 3; Solid -> faces[5][2] = 5;
  Solid -> faces[6][0] = 3; Solid -> faces[6][1] = 2; Solid -> faces[6][2] = 5;
  Solid -> faces[7][0] = 2; Solid -> faces[7][1] = 1; Solid -> faces[7][2] = 5;

#elif DEGENERATE_RECT_PRISM /*Rectangular prism with more faces/edges than necessary*/
  
  Solid -> points[0][0] = 0; Solid -> points[0][1] = -1; Solid -> points[0][2] = 0;
  Solid -> points[1][0] = 1; Solid -> points[1][1] = -1; Solid -> points[1][2] = 0;
  Solid -> points[2][0] = 2; Solid -> points[2][1] = -1; Solid -> points[2][2] = 0;
  Solid -> points[3][0] = 2; Solid -> points[3][1] = 0; Solid -> points[3][2] = 0;
  Solid -> points[4][0] = 0; Solid -> points[4][1] = 0; Solid -> points[4][2] = 0;
  Solid -> points[5][0] = 0; Solid -> points[5][1] = -1; Solid -> points[5][2] = 1;
  Solid -> points[6][0] = 1; Solid -> points[6][1] = -1; Solid -> points[6][2] = 1;
  Solid -> points[7][0] = 2; Solid -> points[7][1] = -1; Solid -> points[7][2] = 1; 
  Solid -> points[8][0] = 2; Solid -> points[8][1] = 0; Solid -> points[8][2] = 1;
  Solid -> points[9][0] = 1; Solid -> points[9][1] = 0; Solid -> points[9][2] = 1;
  Solid -> points[10][0] = 0; Solid -> points[10][1] = 0; Solid -> points[10][2] = 1;
  Solid -> points[11][0] = 0; Solid -> points[11][1] = -0.5; Solid -> points[11][2] = 1;
  Solid -> points[12][0] = 1; Solid -> points[12][1] = -0.5; Solid -> points[12][2] = 1;

  for(i = 0; i < Solid -> num_faces; i++)
    Solid -> face_size[i] = 4;

  for(i = 0; i < Solid -> num_faces; i++) {
    Solid -> faces[i] = malloc(sizeof(int)*(Solid -> face_size[i]));
    if(Solid -> faces[i] == NULL)
      return(ERROR);
  }

  Solid -> faces[0][0] = 10; Solid -> faces[0][1] = 4; Solid -> faces[0][2] = 0; Solid -> faces[0][3] = 5;
  Solid -> faces[1][0] = 0; Solid -> faces[1][1] = 1; Solid -> faces[1][2] = 6; Solid -> faces[1][3] = 5;
  Solid -> faces[2][0] = 1; Solid -> faces[2][1] = 2; Solid -> faces[2][2] = 7; Solid -> faces[2][3] = 6;
  Solid -> faces[3][0] = 8; Solid -> faces[3][1] = 7; Solid -> faces[3][2] = 2; Solid -> faces[3][3] = 3;
  Solid -> faces[4][0] = 3; Solid -> faces[4][1] = 4; Solid -> faces[4][2] = 10; Solid -> faces[4][3] = 8;
  Solid -> faces[5][0] = 6; Solid -> faces[5][1] = 7; Solid -> faces[5][2] = 8; Solid -> faces[5][3] = 9;
  Solid -> faces[6][0] = 11; Solid -> faces[6][1] = 12; Solid -> faces[6][2] = 9; Solid -> faces[6][3] = 10;
  Solid -> faces[7][0] = 5; Solid -> faces[7][1] = 6; Solid -> faces[7][2] = 12; Solid -> faces[7][3] = 11;
  Solid -> faces[8][0] = 0; Solid -> faces[8][1] = 4; Solid -> faces[8][2] = 3; Solid -> faces[8][3] = 2;
  
#endif
  

  printf("So the faces and points of the polyhedron in order are ...\n");
  for(i = 0; i < Solid -> num_faces; i++) {
    for(j = 0; j < Solid -> face_size[i]; j++) {
      printf("(%g,%g,%g) ", Solid -> points[Solid->faces[i][j]][0], Solid -> points[Solid->faces[i][j]][1], 
	     Solid -> points[Solid->faces[i][j]][2]); /*Print a point per face*/
    }
    printf("\n");
  }
}

void traverse_faces(Body2D_tnode tnode)
{
  int i;

  if(tnode -> parent != NULL) {
    printf("Face %d, Parent %d, level %d, bound [%g,%g,%g,%g], num_points %d, ",tnode->poly_num,tnode->parent->poly_num,tnode->level,
	   tnode->bounding_box[0],tnode->bounding_box[1],tnode->bounding_box[2],tnode->bounding_box[3],tnode->Poly->num_points);
    if(tnode == tnode -> parent -> left)
      printf("left, ");
    else if(tnode == tnode -> parent -> right)
      printf("right, ");
    printf("Lwr[%g,%g,%g,%g], Upr[%g,%g,%g,%g]\n",tnode->lower_bound[0],tnode->lower_bound[1],tnode->lower_bound[2],tnode->lower_bound[3],
	   tnode->upper_bound[0],tnode->upper_bound[1],tnode->upper_bound[2],tnode->upper_bound[3]);
    printf("Stored points are - \n");
    for(i = 0; i < tnode->Poly->num_points; i++)
      printf("(%g,%g,%g) ",tnode->Poly->points[i][0],tnode->Poly->points[i][1],tnode->Poly->points[i][2]);
    printf("\n");
    printf("Edges (in order) are - \n");
    for(i = 0; i < tnode->Poly->num_points; i++) {
      printf("(%g,%g,%g) -> (%g,%g,%g)\n", tnode->Poly->points[tnode->Poly->edges[i][0]][0],
	     tnode->Poly->points[tnode->Poly->edges[i][0]][1],tnode->Poly->points[tnode->Poly->edges[i][0]][2],
	     tnode->Poly->points[tnode->Poly->edges[i][1]][0],tnode->Poly->points[tnode->Poly->edges[i][1]][1],
	     tnode->Poly->points[tnode->Poly->edges[i][1]][2]);       
    }
  }
  else {
    printf("Face %d, level %d, bound [%g,%g,%g,%g], num_points %d, ",tnode->poly_num,tnode->level,tnode->bounding_box[0],tnode->bounding_box[1],
	   tnode->bounding_box[2],tnode->bounding_box[3], tnode->Poly->num_points);
    printf("Lwr[%g,%g,%g,%g], Upr[%g,%g,%g,%g]\n",tnode->lower_bound[0],tnode->lower_bound[1],tnode->lower_bound[2],tnode->lower_bound[3],
	   tnode->upper_bound[0],tnode->upper_bound[1],tnode->upper_bound[2],tnode->upper_bound[3]);
    printf("Stored points are - \n");
    for(i = 0; i < tnode->Poly->num_points; i++)
      printf("(%g,%g,%g) ",tnode->Poly->points[i][0],tnode->Poly->points[i][1],tnode->Poly->points[i][2]);
    printf("\n");
    printf("Edges (in order) are - \n");
    for(i = 0; i < tnode->Poly->num_points; i++) {
      printf("(%g,%g,%g) -> (%g,%g,%g)\n", tnode->Poly->points[tnode->Poly->edges[i][0]][0],
	     tnode->Poly->points[tnode->Poly->edges[i][0]][1],tnode->Poly->points[tnode->Poly->edges[i][0]][2],
	     tnode->Poly->points[tnode->Poly->edges[i][1]][0],tnode->Poly->points[tnode->Poly->edges[i][1]][1],
	     tnode->Poly->points[tnode->Poly->edges[i][1]][2]);       
    }
  }

  if(tnode -> left != NULL)
    traverse_faces(tnode->left);
  if(tnode -> right != NULL)
    traverse_faces(tnode->right);
}
