/** \file sspmisc.c
 * \ingroup plot
 * \brief Miscellaneous functions for Simple Scientific Plotting.
 *
 * \author PA Jacobs
 */

/*
 * ------------------
 * Global Definitions
 * ------------------
 */

#include "../../util/source/compiler.h"
#include "../../geometry/source/geom.h"
#include "ssp.h"

#include <stdio.h>
#include <math.h>

#if (STDLIBH)
#  include <stdlib.h>
#  include <stdarg.h>
#endif

/*-------------------------------------------------------------*/

/*
 * Purpose...
 * -------
 * Allocate memory for the data vectors.
 *
 * Input...
 * -----
 * *V    : pointer to the data_array structure
 * idim  : size of the array (0...idim-1)
 *
 */
int ssp_AllocateVector (struct ssp_data_vector *V,
			int idim)
{
if ((V->x = (double *) malloc(idim * sizeof(double))) == NULL)
   goto Fail;
if ((V->value = (double *) malloc(idim * sizeof(double))) == NULL)
   goto Fail;

V->imin = 0;
V->imax = idim - 1;

/* Success: */
   return (0);

Fail:
   printf ("ssp_AllocateVector: Memory allocation failed");
   exit (1);
   return (1);
}

int ssp_FreeVector (struct ssp_data_vector *V) {
    free(V->value);
    free(V->x);
    return 0;
}

/*-----------------------------------------------------------------*/

/*
 * Purpose...
 * -------
 * Allocate memory for the data arrays used in the contouring
 * routines.
 *
 * Input...
 * -----
 * *F    : pointer to the data_array structure
 * ixdim : size of the array in the x-index direction (0...ixdim-1)
 * iydim : size of the array in the y-index direction (0...iydim-1)
 *
 * Revision...
 * --------
 * 12-Jan-94 : added u,v components to the ssp_data_array
 *
 */
int ssp_AllocateArray (struct ssp_data_array *F,
		       int ixdim, int iydim)
{
int ix;


if ((F->x = (double **) malloc(ixdim * sizeof(double *))) == NULL)
   goto Fail;
if ((F->y = (double **) malloc(ixdim * sizeof(double *))) == NULL)
   goto Fail;
if ((F->u = (double **) malloc(ixdim * sizeof(double *))) == NULL)
   goto Fail;
if ((F->v = (double **) malloc(ixdim * sizeof(double *))) == NULL)
   goto Fail;
if ((F->value = (double **) malloc(ixdim * sizeof(double *))) == NULL)
   goto Fail;

for (ix = 0; ix < ixdim; ++ix)
   {
   if ((F->x[ix] = (double *) malloc(iydim * sizeof(double))) == NULL)
      goto Fail;
   if ((F->y[ix] = (double *) malloc(iydim * sizeof(double))) == NULL)
      goto Fail;
   if ((F->u[ix] = (double *) malloc(iydim * sizeof(double))) == NULL)
      goto Fail;
   if ((F->v[ix] = (double *) malloc(iydim * sizeof(double))) == NULL)
      goto Fail;
   if ((F->value[ix] = (double *) malloc(iydim * sizeof(double))) == NULL)
      goto Fail;
   }

F->ixmin = 0;
F->ixmax = ixdim - 1;

F->iymin = 0;
F->iymax = iydim - 1;


/* Success: */
   return (0);

Fail:
   printf ("ssp_AllocateArray: Memory allocation failed");
   exit (1);

return (0);
}

int ssp_FreeArray(struct ssp_data_array *F) {
    int ix;
    for (ix = 0; ix <= F->ixmax; ++ix) {
	free(F->value[ix]);
	free(F->v[ix]);
	free(F->u[ix]);
	free(F->y[ix]);
	free(F->x[ix]);
    }
    free(F->value);
    free(F->v);
    free(F->u);
    free(F->y);
    free(F->x);
    return 0;
}

/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_VectorExtremes (struct ssp_data_vector *V)

#else

int ssp_VectorExtremes (V)
struct ssp_data_vector *V;

#endif

/*
 * Purpose...
 * -------
 * Find the extreme x, f values in the vector.
 *
 * Input...
 * -----
 * V   : pointer to the vector data-structure (see ssp.h)
 *
 */

{  /* begin ssp_VectorExtremes() */
int    i;
double x, v;

/*
 * Set up some outrageous values.
 */
V->xmin = 1.0e12;
V->xmax = -1.0e12;
V->vmin = 1.0e12;
V->vmax = -1.0e12;

/*
 * Now, search the array.
 */
for (i = V->imin; i <= V->imax; ++i)
   {
   x = V->x[i];
   v = V->value[i];

   if (x < V->xmin) V->xmin = x;
   if (x > V->xmax) V->xmax = x;

   if (v < V->vmin) V->vmin = v;
   if (v > V->vmax) V->vmax = v;
   }

return (0);
}  /* end of ssp_VectorExtremes() */
/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_ArrayExtremes (struct ssp_data_array *F)

#else

int ssp_ArrayExtremes (F)
struct ssp_data_array *F;

#endif

/*
 * Purpose...
 * -------
 * Find the extreme x, y, f values in the array.
 *
 * Input...
 * -----
 * F   : pointer to the array structure (see ssp.h)
 *
 */

{  /* begin ssp_ArrayExtremes() */
int    ix, iy;
double x, y, v;

/*
 * Set up some outrageous values.
 */
F->xmin = 1.0e12;
F->xmax = -1.0e12;
F->ymin = 1.0e12;
F->ymax = -1.0e12;
F->vmin = 1.0e12;
F->vmax = -1.0e12;

/*
 * Now, search the array.
 */
for (ix = F->ixmin; ix <= F->ixmax; ++ix)
   for (iy = F->iymin; iy <= F->iymax; ++iy)
      {
      x = F->x[ix][iy];
      y = F->y[ix][iy];
      v = F->value[ix][iy];

      if (x < F->xmin) F->xmin = x;
      if (x > F->xmax) F->xmax = x;

      if (y < F->ymin) F->ymin = y;
      if (y > F->ymax) F->ymax = y;

      if (v < F->vmin) F->vmin = v;
      if (v > F->vmax) F->vmax = v;
      }

return (0);
}  /* end of ssp_ArrayExtremes() */
/*-----------------------------------------------------------------*/
