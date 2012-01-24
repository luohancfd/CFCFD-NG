#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gio_kernel.h"
#include "gio_adts.h"

#define SQRT(A) (pow((A), 0.5))

#define TRIANGULAR_PRISM 0
#define DEGENERATE_TETRAHEDRON 1

void get_point(double []);
void traverse_faces(Body2D_tnode);
short int make_polyh(Body);

main()
{
  Body Solid;
  double point[3];
  short int i;

  Solid = malloc(sizeof(struct body));
  if(Solid == NULL)
    return(ERROR);

  /*Define polyhedron*/
  get_point(point);
  make_polyh(Solid);

  /*Build ADT*/
  if(build_face_ADT(Solid) == ERROR)
    return(ERROR);

  traverse_faces(Solid -> Faces);
  printf("\n");

  i = polyhio1(point, Solid);
  if(i == IN_POLY)
    printf("\nINSIDE\n");
  else if(i == OUT_POLY)
    printf("\nOUTSIDE\n");
}

void get_point(double point[]) /*Define points here for the moment*/
{
  point[0] = 2 ;
  point[1] = SQRT(3.0)/2  ;
  point[2] = SQRT(3.0) + 2  ;

  printf("The point I'm queryin' is (%g,%g,%g)\n",point[0],point[1],point[2]);
}

short int make_polyh(Body Solid)
{
  int i, j;

#if TRIANGULAR_PRISM
  Solid -> num_points = 6;
  Solid -> num_faces = 5;

#elif DEGENERATE_TETRAHEDRON
  Solid -> num_points = 5;
  Solid -> num_faces = 5;
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

#if TRIANGULAR_PRISM
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

#endif
  

  printf("So the faces and points of the polyhedron in order are ...\n");
  for(i = 0; i < Solid -> num_faces; i++) {
    for(j = 0; j < Solid -> face_size[i]; j++) {
      printf("(%g,%g,%g) ", Solid->points[Solid->faces[i][j]][0], Solid -> points[Solid->faces[i][j]][1], 
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
