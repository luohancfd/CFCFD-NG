/** \file bspatch.c
 * \brief B-spline patches
 * \ingroup nm
 *
 *-----------------------------------------------------------*
 *       bs_init_arrays   ->  reads control net data         *
 *       bs_calc          ->  interpolates parametric point  *
 *-----------------------------------------------------------*
 *
 * Purpose...
 *
 * Read in the control points for 4 B-spline patches and generate
 * an interpolated point within the volume enclosed by the 4 patches.
 *
 * \author C. Craddock
 *
 * \date September 1993
 *
 * \version 29-Sep-93 : Adapted for the "sm_3d" Space-marching flow solver (PJ)
 * \version 30-Sep-94 : Make the number of knots independent for each
 *             bspline surface
 * \version 09-Mar-94 : Increase the size of the workspace arrays 
 *
 */
#include <stdio.h>
#include <math.h>

#ifndef COMPILER_H
#  include "../../util/source/compiler.h"
#endif

#if (STDLIBH)
#  include <stdlib.h>
#  include <string.h>
#endif

#include "../../sm_3d_plus/source/sm_3d.h"
#include "bspatch.h"

/*---------------------------------------------------------------*/
/* Global Variables */

int     nx1, ny1, nx2, ny2, nx3, ny3, nx4, ny4;

struct  point_3D  **q1, **q2, **q3, **q4;
double  temp1, temp2, temp3;

/*------------------------------------------------------------------*/

/*
 *  function name: bs_init_arrays()
 *
 *  This program reads the contol net data from the appropriate
 *  files and assigns it to structured arrays ready for processing.
 *  The control nets are B-Spline nets that are created by the
 *  program "invsurf".
 *
 *  Written by Chris Craddock
 *  Version 1.0  6/9/93
 *
 *
 */

#if PROTO

int bs_init_arrays ( char base_file_name[32] )

#else

int bs_init_arrays ( base_file_name )
char base_file_name[32];

#endif

{  /* begin bs_init_arrays() */

int i,j;
int failed;  /* flag for failed memory allocation */
char file_name[32];
/*
 *  open the control net files for each patch
 */

FILE *net1, *net2, *net3, *net4;

/*
 * Read control net data
 */

printf ("bs_init_arrays(): reading control net 1\n");

strcpy (file_name, base_file_name); 
strcat (file_name, ".q1");
net1 = NULL;
if ((net1 = fopen (file_name,"r")) == NULL)
   {
   printf ("bs_init_arrays(): Couldn't open q1\n");
   exit(-1);
   }

fscanf (net1, "%d %d", &ny1, &nx1);

failed = 0;
q1 = (struct point_3D **) malloc ( (ny1+1) * sizeof(struct point_3D *) );
if (q1 == NULL) ++failed;
for (i = 0; i <= ny1; ++i)
   {
   q1[i] = (struct point_3D *) malloc ( (nx1+1) * sizeof(struct point_3D) );
   if (q1[i] == NULL) ++failed;
   }

if (failed > 0)
   {
   printf ("bs_init_arrays() : Could not allocate memory.\n");
   exit (-1);
   }

for (j = 0; j <= nx1; j++)
    {
    for (i = 0; i <= ny1; i++)
        {
	fscanf (net1,"%lf %lf %lf",
	   &(q1[i][j].x), &(q1[i][j].y), &(q1[i][j].z) );
        }
    }

fclose(net1);

printf ("bs_init_arrays(): reading control net 2\n");

strcpy (file_name, base_file_name);
strcat (file_name, ".q2");
net2 = NULL;
if ((net2 = fopen (file_name,"r")) == NULL)
   {
   printf ("bs_init_arrays(): Couldn't open q2\n");
   exit(-1);
   }

fscanf (net2,"%d %d",&ny2,&nx2);

failed = 0;
q2 = (struct point_3D **) malloc ( (ny2+1) * sizeof(struct point_3D *) );
if (q2 == NULL) ++failed;
for (i = 0; i <= ny2; ++i)
   {
   q2[i] = (struct point_3D *) malloc ( (nx2+1) * sizeof(struct point_3D) );
   if (q2[i] == NULL) ++failed;
   }

if (failed > 0)
   {
   printf ("bs_init_arrays() : Could not allocate memory.\n");
   exit (-1);
   }

for (j = 0; j <= nx2; j++)
    {
    for (i = 0; i <= ny2; i++)
        {
	fscanf (net2,"%lf %lf %lf",
	   &(q2[i][j].x), &(q2[i][j].y), &(q2[i][j].z) );
        }
    }

fclose(net2);

printf ("bs_init_arrays(): reading control net 3\n");

strcpy (file_name, base_file_name);
strcat (file_name, ".q3");
net3 = NULL;
if ((net3 = fopen (file_name,"r")) == NULL)
   {
   printf ("bs_init_arrays(): Couldn't open q3\n");
   exit(-1);
   }

fscanf (net3,"%d %d",&ny3,&nx3);

failed = 0;
q3 = (struct point_3D **) malloc ( (ny3+1) * sizeof(struct point_3D *) );
if (q3 == NULL) ++failed;
for (i = 0; i <= ny3; ++i)
   {
   q3[i] = (struct point_3D *) malloc ( (nx3+1) * sizeof(struct point_3D) );
   if (q3[i] == NULL) ++failed;
   }

if (failed > 0)
   {
   printf ("bs_init_arrays() : Could not allocate memory.\n");
   exit (-1);
   }

for (j = 0; j <= nx3; j++)
    {
    for (i = 0; i <= ny3; i++)
        {
	fscanf (net3,"%lf %lf %lf",
	   &(q3[i][j].x), &(q3[i][j].y), &(q3[i][j].z) );
        }
    }

fclose(net3);

printf ("bs_init_arrays(): reading control net 4\n");

strcpy (file_name, base_file_name);
strcat (file_name, ".q4");
net4 = NULL;
if ((net4 = fopen (file_name,"r")) == NULL)
   {
   printf ("bs_init_arrays(): Couldn't open q4\n");
   exit(-1);
   }

fscanf (net4,"%d %d",&ny4,&nx4);

failed = 0;
q4 = (struct point_3D **) malloc ( (ny4+1) * sizeof(struct point_3D *) );
if (q4 == NULL) ++failed;
for (i = 0; i <= ny4; ++i)
   {
   q4[i] = (struct point_3D *) malloc ( (nx4+1) * sizeof(struct point_3D) );
   if (q4[i] == NULL) ++failed;
   }

if (failed > 0)
   {
   printf ("bs_init_arrays() : Could not allocate memory.\n");
   exit (-1);
   }

for (j = 0; j <= nx4; j++)
    {
    for (i = 0; i <= ny4; i++)
        {
	fscanf (net4,"%lf %lf %lf",
	   &(q4[i][j].x), &(q4[i][j].y), &(q4[i][j].z) );
        }
    }

fclose(net4);

return 0;
} /* End of function bs_init_arrays */


/*---------------------------------------------------------------------*/

/*
 *  function name: bs_calc(neta,nzeta,xi,eta,zeta,ed1,ed2,ed3,ed4)
 *
 *  Given a parametric position inside the inlet, computes the
 *  three dimensional coordinates of the point. The control nets
 *  are B-Spline nets previously read in by bs_init_arrays()
 *
 *  Written by Chris Craddock
 *  Version 1.0   6/ 9/93  Calculates all points within slice
 *  Version 2.0  29/10/93  Calculates points on edges of slice
 *
 *
 */

#if PROTO

int bs_calc ( int neta, int nzeta,
              double xi, double eta[], double zeta[],
              struct point_3D ed1[],
              struct point_3D ed2[],
              struct point_3D ed3[],
              struct point_3D ed4[])

#else

int bs_calc ( neta, nzeta, xi, eta, zeta, ed1, ed2, ed3, ed4 )
int    neta, nzeta
double xi, eta[], zeta[];
struct point_3D ed1[], ed2[], ed3[], ed4[];

#endif

{  /* begin bs_calc() */

int     v, w;
double  i, etain[NDIM], zetain3[NDIM], zetain4[NDIM];
double  nxend, nyend, nzend3, nzend4;

#if DEBUG_BSPATCH >= 1
   printf ("bs_calc(): Start of function...\n");
   fflush (stdout);
#endif

nxend = nx1 - 2;
nyend = ny1 - 2;
nzend3 = ny3 - 2;
nzend4 = ny4 - 2;


/* 
 * Convert parametric coordinates xi,eta,zeta to relative coord
 */

i = xi   * nxend;

for (v = 0 ; v <= neta ; v++ )
    {
    etain[v] = eta[v] * nyend;
    }

for (w = 0 ; w <= nzeta ; w++ )
    {
    zetain3[w] = zeta[w] * nzend3;
    zetain4[w] = zeta[w] * nzend4;
    }


/*
 * Calculate surface edge points on each patch corresponding
 * to the parametric coordinates
 */

for ( v = 0 ; v <= neta ; v++ )
    bs_calc2( q1, etain[v], i, nx1, ny1, &(ed1[v]) );

for ( v = 0 ; v <= neta ; v++ )
    bs_calc2( q2, etain[v], i, nx2, ny2, &(ed2[v]) );

for ( w = 0 ; w <= nzeta ; w++ )
    bs_calc2( q3, zetain3[w], i, nx3, ny3, &(ed3[w]) );

for ( w = 0 ; w <= nzeta ; w++ )
    bs_calc2( q4, zetain4[w], i, nx4, ny4, &(ed4[w]) );

return (0);
}  /* End of bs_calc */


/*---------------------------------------------------------------*/

/*
 *  The function, bs_calc2 , calculates the point on a surface
 *  patch corresponding to the passed parametric coordinate
 *  u,v
 */

#if PROTO

int bs_calc2 ( struct point_3D **q,
               double u, double v,
               int nx, int ny,
               struct point_3D *pnt )

#else

int bs_calc2 ( q, u, v, nx, ny, pnt )
int    nx, ny;
struct point_3D **q;
double u, v;
struct point_3D *pnt;

#endif

{  /* begin bs_calc2() */

/*
 * Declare internal variables
 */

int    m,n;
double nreal,mreal;
double sumx[NDIM], sumy[NDIM], sumz[NDIM];
double sumxx, sumyy, sumzz;
double fn;

#if DEBUG_BSPATCH >= 1
   printf ("bs_calc2(): Begin function...\n");
   printf ("bs_calc2(): u = %f, v = %f\n", u, v);
   fflush (stdout);
#endif

/*
 * Perform first stage of calculations along x axis
 */

for (m = 0; m <= ny; m++)
    {
    sumx[m] = 0.0;
    sumy[m] = 0.0;
    sumz[m] = 0.0;
    for (n = 0; n <= nx; n++)
        {
        nreal = n;
        fn = bs_nval (nreal, v);
        sumx[m] += q[m][n].x * fn;
        sumy[m] += q[m][n].y * fn;
        sumz[m] += q[m][n].z * fn;
        }
    }

/*
 *  Perform second stage of calculations along y axis
 */
sumxx = 0.0;
sumyy = 0.0;
sumzz = 0.0;

for (m = 0 ; m <= ny ; m++)
    {
    mreal = m;
    fn = bs_nval (mreal, u);
    sumxx += sumx[m] * fn;
    sumyy += sumy[m] * fn;
    sumzz += sumz[m] * fn;
    }

pnt->x = sumxx;
pnt->y = sumyy;
pnt->z = sumzz;

return 0;
}  /* end of bs_calc2() */



/*----------------------------------------------------------------*/

/*
 *  this function calculates the value of a fourth order
 *  b-spline function given x and t, where x is the i,j
 *  ordinate and t is the u,w value.
 */

#if PROTO

double bs_nval ( double x, double t )

#else

double bs_nval (x,t)
double x,t;

#endif

{  /* begin bs_nval() */
double res;


if ( t < (x-3.0) || t >= (x+1.0))
   {res = 0.0;}
else
   {
   if ( t < (x-2.0) )
      {
      temp1 = (t+3.0-x);
      res =  ((1.0/6.0)*temp1*temp1*temp1);
      }
   else
      {
      if ( t < (x-1.0) )
         {
         temp1 = (t-x);
         temp2 = (t+1.0-x);
         temp3 = (t-1.0-x);
         res = ( (2.0/3.0)*temp1*temp1*temp1 - temp2*temp2*temp2 -
               (1.0/6.0)*temp3*temp3*temp3 );
         }
      else
         {
         if ( t < x )
            {
            temp1 = (t-x);
            temp2 = (t-1.0-x);
            res = ( (2.0/3.0)*temp1*temp1*temp1 -
                  (1.0/6.0)*temp2*temp2*temp2 );
            }
         else
            {
            temp1 = (t-1.0-x);
            res = (0.0 - (1.0/6.0)*temp1*temp1*temp1 );
            }
         }
      }
   }

return (res);
}  /* end of bs_nval() */


/*------------------------ end of bspatch.c ------------------------*/
