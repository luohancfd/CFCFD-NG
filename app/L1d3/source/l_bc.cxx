/** \file l_bc.cxx
 * \ingroup l1d2
 * \brief Boundary Conditions for l1d.c
 * 
 * \author PA Jacobs
 */


/*-----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../../lib/util/source/useful.h"
#include "l1d.hh"
#include "l_cell.hh"
#include "l_slug.hh"


int L_bc_left_velocity(GasSlug* A, double v)
{
    /* 
     * Apply a specified velocity to the left boundary.
     */
    int ix;
    LCell *src, *dest;
    double xbndy;

    ix = A->ixmin - 1;
    xbndy = A->Cell[ix].x;
    A->left_ustar = v;

    dest = &(A->Cell[ix]);
    src = &(A->Cell[ix + 1]);
    L_copy_cell_data(src, dest, 0);
    dest->u = v;
    dest->xmid = xbndy - (src->xmid - xbndy);

    dest = &(A->Cell[ix - 1]);
    src = &(A->Cell[ix + 2]);
    L_copy_cell_data(src, dest, 0);
    dest->u = v;
    dest->xmid = xbndy - (src->xmid - xbndy);

    return 0;
}   /* end function L_bc_left_velocity */


int L_bc_left_reflect(GasSlug* A)
{
    /*
     * Apply a reflective bc to the left boundary.
     */
    int ix;
    LCell *src, *dest;
    double xbndy;

    ix = A->ixmin - 1;
    xbndy = A->Cell[ix].x;
    A->left_ustar = 0.0;

    dest = &(A->Cell[ix]);
    src = &(A->Cell[ix + 1]);
    L_copy_cell_data(src, dest, 0);
    dest->u *= -1.0;
    dest->xmid = xbndy - (src->xmid - xbndy);

    dest = &(A->Cell[ix - 1]);
    src = &(A->Cell[ix + 2]);
    L_copy_cell_data(src, dest, 0);
    dest->u *= -1.0;
    dest->xmid = xbndy - (src->xmid - xbndy);

    return 0;
}   /* end function L_bc_left_reflect */


int L_bc_left_free(GasSlug *A)
{
    /*
     * Apply a free/supersonic outflow bc to the left boundary.
     */
    int ix;
    LCell *src, *dest;
    double xbndy;

    ix = A->ixmin - 1;
    xbndy = A->Cell[ix].x;
    /*
     * We should not need this velocity, but we'll copy it anyway.
     */
    A->left_ustar = A->Cell[ix + 1].u;

    dest = &(A->Cell[ix]);
    src = &(A->Cell[ix + 1]);
    L_copy_cell_data(src, dest, 0);
    dest->xmid = xbndy - (src->xmid - xbndy);

    dest = &(A->Cell[ix - 1]);
    src = &(A->Cell[ix + 1]);
    L_copy_cell_data(src, dest, 0);
    dest->xmid = xbndy - (A->Cell[ix + 2].xmid - xbndy);

    return 0;
}   /* end function L_bc_left_free */


int L_bc_right_velocity(GasSlug* A, double v)
{
    /*
     * Apply a specified velocity to the right boundary.
     */
    int ix;
    LCell *src, *dest;
    double xbndy;

    ix = A->ixmax + 1;
    xbndy = A->Cell[ix - 1].x;
    A->right_ustar = v;

    dest = &(A->Cell[ix]);
    src = &(A->Cell[ix - 1]);
    L_copy_cell_data(src, dest, 0);
    dest->u = v;
    dest->xmid = xbndy + (xbndy - src->xmid);

    dest = &(A->Cell[ix + 1]);
    src = &(A->Cell[ix - 2]);
    L_copy_cell_data(src, dest, 0);
    dest->u = v;
    dest->xmid = xbndy + (xbndy - src->xmid);

    return 0;
}   /* end function L_bc_right_velocity */



int L_bc_right_reflect(GasSlug* A)
{
    /*
     * Apply a reflective bc to the right boundary.
     */
    int ix;
    LCell *src, *dest;
    double xbndy;

    ix = A->ixmax + 1;
    xbndy = A->Cell[ix - 1].x;
    A->right_ustar = 0.0;

    dest = &(A->Cell[ix]);
    src = &(A->Cell[ix - 1]);
    L_copy_cell_data(src, dest, 0);
    dest->u *= -1.0;
    dest->xmid = xbndy + (xbndy - src->xmid);

    dest = &(A->Cell[ix + 1]);
    src = &(A->Cell[ix - 2]);
    L_copy_cell_data(src, dest, 0);
    dest->u *= -1.0;
    dest->xmid = xbndy + (xbndy - src->xmid);

    return 0;
}   /* end function L_bc_right_reflect */


int L_bc_right_free(GasSlug* A)
{
    /* 
     * Apply a free/supersonic-outflow bc to the right boundary.
     */
    int ix;
    LCell *src, *dest;
    double xbndy;

    ix = A->ixmax + 1;
    xbndy = A->Cell[ix - 1].x;
    /*
     * We should not need this velocity, but we'll copy it anyway. 
     */
    A->right_ustar = A->Cell[ix - 1].u;

    dest = &(A->Cell[ix]);
    src = &(A->Cell[ix - 1]);
    L_copy_cell_data(src, dest, 0);
    dest->xmid = xbndy + (xbndy - src->xmid);

    dest = &(A->Cell[ix + 1]);
    src = &(A->Cell[ix - 2]);
    L_copy_cell_data(src, dest, 0);
    dest->xmid = xbndy + (xbndy - src->xmid);

    return 0;
}   /* end function L_bc_right_free */


int L_exchange_bc_data(GasSlug* A, GasSlug* B)
{
    /*
     * Copy the end data from the left slug (A) to the  right slug (B).
     *       ia-1 ia    |  ia+1 ia+2
     *        V    V    |   ^    ^
     *       ib-2 ib-1  |  ib   ib+1
     */
    int ia, ib;
    LCell *src, *dest;

    ia = A->ixmax;
    ib = B->ixmin;

    A->right_ustar = B->Cell[ib].u;
    B->left_ustar = A->Cell[ia].u;

    dest = &(A->Cell[ia + 1]);
    src = &(B->Cell[ib]);
    L_copy_cell_data(src, dest, 0);

    dest = &(A->Cell[ia + 2]);
    src = &(B->Cell[ib + 1]);
    L_copy_cell_data(src, dest, 0);

    dest = &(B->Cell[ib - 1]);
    src = &(A->Cell[ia]);
    L_copy_cell_data(src, dest, 0);

    dest = &(B->Cell[ib - 2]);
    src = &(A->Cell[ia - 1]);
    L_copy_cell_data(src, dest, 0);

    return 0;
}   /* end function L_exchange_bc_data() */


int L_blend_slug_ends(GasSlug* A, int endA, 
		      GasSlug* B, int endB,
		      double dxb )
{
    /*
     * Purpose...
     * -------
     * Blend the ends of the specified gas slugs over the distance dxb.
     * Choose the type of blending below.
     *
     * Input...
     * -----
     * *A *B     : pointer to the gas slug data structures.
     * endA, endB: integers to indicate which end of the gas slugs we blend
     * dxb       : distance over which we will blend the flow properties
     *
     * Returns...
     * -------
     * 0 if nothing drastic happens.
     */
    int ix, iA, iB, blend_type;
    LCell *c, *cA, *cB;
    double x, x0, xA, xB, half_dxb, alpha;

    blend_type = BLEND_PUT;
    /* blend_type = BLEND_RHOUE; */

#   if DEBUG >= 2
    printf( "\nBlend end of slugs over distance %g with type %d...\n", 
	    dxb, blend_type );
#   endif

    if ( dxb <= 0.0 ) return 0;

    /*
     * Search from the ends of the slugs to find the indices over which
     * we will blend.
     * If dxb/2 is larger than the whole slug, we'll just end up blending
     * over the whole slug; this should be safe enough.
     */
    xA = 0.0;
    xB = 0.0;
    /* gcc could not determine if the following code initialized xA, xB */
    half_dxb = 0.5 * dxb;
    if ( endA == LEFT ) {
        iA = A->ixmin;
        x0 = A->Cell[A->ixmin - 1].x;
        for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
            c = &(A->Cell[ix]);
            xA = c->xmid;
	    iA = ix;
	    if ( fabs(xA - x0) > half_dxb ) break;
        }   /* end for */
    } else {
        /* Assume right-hand end. */
        iA = A->ixmax;
        x0 = A->Cell[A->ixmax].x;
        for (ix = A->ixmax; ix >= A->ixmin; --ix) {
            c = &(A->Cell[ix]);
            xA = c->xmid;
	    iA = ix;
	    if ( fabs(xA - x0) > half_dxb ) break;
        }   /* end for */
    }   /* end if */

    if ( endB == LEFT ) {
        iB = B->ixmin;
        x0 = B->Cell[B->ixmin - 1].x;
        for (ix = B->ixmin; ix <= B->ixmax; ++ix) {
            c = &(B->Cell[ix]);
            xB = c->xmid;
	    iB = ix;
	    if ( fabs(xB - x0) > half_dxb ) break;
        }   /* end for */
    } else {
        /* Assume right-hand end. */
        iB = B->ixmax;
        x0 = B->Cell[B->ixmax].x;
        for (ix = B->ixmax; ix >= B->ixmin; --ix) {
            c = &(B->Cell[ix]);
            xB = c->xmid;
	    iB = ix;
	    if ( fabs(xB - x0) > half_dxb ) break;
        }   /* end for */
    }   /* end if */

#   if DEBUG >= 2
    printf( "Blend slugs between cells %d and %d over distance %g\n", 
	    iA, iB, fabs(xB - xA) );
#   endif

    cA = &( A->Cell[iA] );
    cB = &( B->Cell[iB] );
    if ( endA == LEFT ) {
        for (ix = A->ixmin; ix < iA; ++ix) {
            c = &(A->Cell[ix]);
            x = c->xmid;
            alpha = (x - xA) / (xB - xA);
            L_blend_cells( cA, cB, c, alpha, blend_type );
        }   /* end for */
    } else {
        for (ix = A->ixmax; ix > iA; --ix) {
            c = &(A->Cell[ix]);
            x = c->xmid;
            alpha = (x - xA) / (xB - xA);
            L_blend_cells( cA, cB, c, alpha, blend_type );
        }   /* end for */
    }

    if ( endB == LEFT ) {
        for (ix = B->ixmin; ix < iB; ++ix) {
            c = &(B->Cell[ix]);
            x = c->xmid;
            alpha = (x - xA) / (xB - xA);
            L_blend_cells( cA, cB, c, alpha, blend_type );
        }   /* end for */
    } else {
        for (ix = B->ixmax; ix > iB; --ix) {
            c = &(B->Cell[ix]);
            x = c->xmid;
            alpha = (x - xA) / (xB - xA);
            L_blend_cells( cA, cB, c, alpha, blend_type );
        }   /* end for */
    }

    return 0;
}   /* end function L_blend_slug_ends() */
