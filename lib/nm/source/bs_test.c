/** \file bs_test.c
 * \ingroup nm
 * \brief Test code for bs_* series of codes
 */

#include "../../util/source/compiler.h"
#include "../../plot/source/geom.h"
#include "bspatch.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INC 0.1     /* Deteremines the number of grid points */

#define NDIM 61
struct point_3D ed1[NDIM], ed2[NDIM], ed3[NDIM], ed4[NDIM];

main ()
{
int    neta,nzeta;
double xi,eta[NDIM],zeta[NDIM],v,w;

/*
 *  Open file to contatian slice data
 */

FILE *grid;
grid = fopen("inlet.xyz","w");

bs_init_arrays ();

/*
 *  Get slice position
 */

printf ("Please enter xi value indicating slice position :\n");
scanf("%lf", &xi);

/*
 *  Determine parametric points
 */

neta = -1;
nzeta = -1;
v = 0.0;

for ( v = 0 ; v <= 1.00001  ; v = v + INC )
    {
    neta = neta + 1;
    eta[neta] = v ;
    }

for ( w = 0 ; w <= 1.00001  ; w = w + INC )
    {
    nzeta = nzeta + 1;
    zeta[nzeta] = w;
    }

/*
 *  Call function bs_calc
 */

bs_calc (neta, nzeta, xi, eta, zeta, ed1, ed2, ed3, ed4);

/*
 *  Write data to file called inlet.xyz
 */

neta = -1;
fprintf (grid, "Side 1 \n");

for ( v = 0 ; v <= 1.00001  ; v = v + INC)
    {
    neta = neta + 1;
    fprintf (grid, "%lf %lf %lf \n ",ed1[neta].x,ed1[neta].y,ed1[neta].z);
    }

neta = -1;
fprintf (grid, "Side 2 \n");

for ( v = 0 ; v <= 1.00001 ; v = v + INC)
    {
    neta = neta + 1;
    fprintf (grid, "%lf %lf %lf \n ",ed2[neta].x,ed2[neta].y,ed2[neta].z);
    }

nzeta = -1;
fprintf (grid, "Side 3 \n");

for ( w = 0 ; w <= 1.00001  ; w = w + INC)
    {
    nzeta = nzeta + 1;
    fprintf (grid, "%lf %lf %lf \n ",ed3[nzeta].x,ed3[nzeta].y,ed3[nzeta].z);
    }

nzeta = -1;
fprintf (grid, "Side 4 \n");

for ( w = 0 ; w <= 1.00001  ; w = w + INC)
    {
    nzeta = nzeta + 1;
    fprintf (grid, "%lf %lf %lf \n ",ed4[nzeta].x,ed4[nzeta].y,ed4[nzeta].z);
    }

fclose(grid);

return 0;

/*  End of main */
}

