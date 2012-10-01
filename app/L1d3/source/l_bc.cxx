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
    // Apply a specified velocity to the left boundary.
    int ix = A->ixmin - 1;
    double xbndy = A->Cell[ix].x;
    A->left_ustar = v;
    LCell& dest =A->Cell[ix];
    LCell& src = A->Cell[ix+1];
    dest.copy_data_from(src, 0);
    dest.u = v;
    dest.xmid = xbndy - (src.xmid - xbndy);
    LCell& dest2 = A->Cell[ix-1];
    LCell& src2 = A->Cell[ix+2];
    dest2.copy_data_from(src2, 0);
    dest2.u = v;
    dest2.xmid = xbndy - (src2.xmid - xbndy);
    return SUCCESS;
}   /* end function L_bc_left_velocity */


int L_bc_left_reflect(GasSlug* A)
{
    // Apply a reflective bc to the left boundary.
    int ix = A->ixmin - 1;
    double xbndy = A->Cell[ix].x;
    A->left_ustar = 0.0;
    LCell& dest = A->Cell[ix];
    LCell& src = A->Cell[ix+1];
    dest.copy_data_from(src, 0);
    dest.u *= -1.0;
    dest.xmid = xbndy - (src.xmid - xbndy);
    LCell& dest2 = A->Cell[ix-1];
    LCell& src2 = A->Cell[ix+2];
    dest2.copy_data_from(src2, 0);
    dest2.u *= -1.0;
    dest2.xmid = xbndy - (src2.xmid - xbndy);
    return SUCCESS;
}   /* end function L_bc_left_reflect */


int L_bc_left_free(GasSlug *A)
{
    // Apply a free/supersonic outflow bc to the left boundary.
    int ix = A->ixmin - 1;
    double xbndy = A->Cell[ix].x;
    A->left_ustar = A->Cell[ix+1].u;
    LCell& dest = A->Cell[ix];
    LCell& src = A->Cell[ix+1];
    dest.copy_data_from(src, 0);
    dest.xmid = xbndy - (src.xmid - xbndy);
    LCell& dest2 = A->Cell[ix-1];
    LCell& src2 = A->Cell[ix+2];
    dest2.copy_data_from(src, 0); // same as for first ghost cell
    dest2.xmid = xbndy - (src2.xmid - xbndy);
    return SUCCESS;
}   /* end function L_bc_left_free */


int L_bc_right_velocity(GasSlug* A, double v)
{
    // Apply a specified velocity to the right boundary.
    int ix = A->ixmax + 1;
    double xbndy = A->Cell[ix - 1].x;
    A->right_ustar = v;
    LCell& dest = A->Cell[ix];
    LCell& src = A->Cell[ix-1];
    dest.copy_data_from(src, 0);
    dest.u = v;
    dest.xmid = xbndy + (xbndy - src.xmid);
    LCell& dest2 = A->Cell[ix+1];
    LCell& src2 = A->Cell[ix-2];
    dest2.copy_data_from(src2, 0);
    dest2.u = v;
    dest2.xmid = xbndy + (xbndy - src2.xmid);
    return SUCCESS;
}   /* end function L_bc_right_velocity */



int L_bc_right_reflect(GasSlug* A)
{
    // Apply a reflective bc to the right boundary.
    int ix = A->ixmax + 1;
    double xbndy = A->Cell[ix - 1].x;
    A->right_ustar = 0.0;
    LCell& dest = A->Cell[ix];
    LCell& src = A->Cell[ix-1];
    dest.copy_data_from(src, 0);
    dest.u *= -1.0;
    dest.xmid = xbndy + (xbndy - src.xmid);
    LCell& dest2 = A->Cell[ix+1];
    LCell& src2 = A->Cell[ix-2];
    dest2.copy_data_from(src2, 0);
    dest2.u *= -1.0;
    dest2.xmid = xbndy + (xbndy - src2.xmid);
    return SUCCESS;
}   /* end function L_bc_right_reflect */


int L_bc_right_free(GasSlug* A)
{
    // Apply a free/supersonic-outflow bc to the right boundary.
    int ix = A->ixmax + 1;
    double xbndy = A->Cell[ix - 1].x;
    A->right_ustar = A->Cell[ix - 1].u;
    LCell& dest = A->Cell[ix];
    LCell& src = A->Cell[ix-1];
    dest.copy_data_from(src, 0);
    dest.xmid = xbndy + (xbndy - src.xmid);
    LCell& dest2 = A->Cell[ix+1];
    LCell& src2 = A->Cell[ix-2];
    dest2.copy_data_from(src, 0); // same as for first ghost cell
    dest2.xmid = xbndy + (xbndy - src2.xmid);
    return SUCCESS;
}   /* end function L_bc_right_free */


int L_exchange_bc_data(GasSlug* A, GasSlug* B)
{
    /*
     * Copy the end data from the left slug (A) to the  right slug (B).
     *       ia-1 ia    |  ia+1 ia+2
     *        V    V    |   ^    ^
     *       ib-2 ib-1  |  ib   ib+1
     */
    int ia = A->ixmax;
    int ib = B->ixmin;
    A->right_ustar = B->Cell[ib].u;
    B->left_ustar = A->Cell[ia].u;
    A->Cell[ia+1].copy_data_from(B->Cell[ib], 0);
    A->Cell[ia+2].copy_data_from(B->Cell[ib+1], 0);
    B->Cell[ib-1].copy_data_from(A->Cell[ia], 0);
    B->Cell[ib-2].copy_data_from(A->Cell[ia-1], 0);
    return SUCCESS;
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
     * *A *B     : pointer to the GasSlugs.
     * endA, endB: integers to indicate which end of the gas slugs we blend
     * dxb       : distance over which we will blend the flow properties
     *
     * Returns...
     * -------
     * 0 if nothing drastic happens.
     */
    int ix, iA, iB, blend_type;
    double x, x0, xA, xB, half_dxb, alpha;

    blend_type = BLEND_PUT;
    /* blend_type = BLEND_RHOUE; */
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
        x0 = A->Cell[A->ixmin-1].x;
        for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
            xA = A->Cell[ix].xmid;
	    iA = ix;
	    if ( fabs(xA - x0) > half_dxb ) break;
        }   /* end for */
    } else {
        /* Assume right-hand end. */
        iA = A->ixmax;
        x0 = A->Cell[A->ixmax].x;
        for (ix = A->ixmax; ix >= A->ixmin; --ix) {
            xA = A->Cell[ix].xmid;
	    iA = ix;
	    if ( fabs(xA - x0) > half_dxb ) break;
        }   /* end for */
    }   /* end if */

    if ( endB == LEFT ) {
        iB = B->ixmin;
        x0 = B->Cell[B->ixmin - 1].x;
        for (ix = B->ixmin; ix <= B->ixmax; ++ix) {
            xB = B->Cell[ix].xmid;
	    iB = ix;
	    if ( fabs(xB - x0) > half_dxb ) break;
        }   /* end for */
    } else {
        /* Assume right-hand end. */
        iB = B->ixmax;
        x0 = B->Cell[B->ixmax].x;
        for (ix = B->ixmax; ix >= B->ixmin; --ix) {
            xB = B->Cell[ix].xmid;
	    iB = ix;
	    if ( fabs(xB - x0) > half_dxb ) break;
        }   /* end for */
    }   /* end if */

    LCell& cA = A->Cell[iA];
    LCell& cB = B->Cell[iB];
    if ( endA == LEFT ) {
        for (ix = A->ixmin; ix < iA; ++ix) {
            LCell& c = A->Cell[ix];
            x = c.xmid;
            alpha = (x - xA) / (xB - xA);
            L_blend_cells(cA, cB, c, alpha, blend_type);
        }   /* end for */
    } else {
        for (ix = A->ixmax; ix > iA; --ix) {
            LCell& c = A->Cell[ix];
            x = c.xmid;
            alpha = (x - xA) / (xB - xA);
            L_blend_cells(cA, cB, c, alpha, blend_type);
        }   /* end for */
    }

    if ( endB == LEFT ) {
        for (ix = B->ixmin; ix < iB; ++ix) {
            LCell& c = B->Cell[ix];
            x = c.xmid;
            alpha = (x - xA) / (xB - xA);
            L_blend_cells(cA, cB, c, alpha, blend_type);
        }   /* end for */
    } else {
        for (ix = B->ixmax; ix > iB; --ix) {
            LCell& c = B->Cell[ix];
            x = c.xmid;
            alpha = (x - xA) / (xB - xA);
            L_blend_cells(cA, cB, c, alpha, blend_type);
        }   /* end for */
    }
    return SUCCESS;
}   /* end function L_blend_slug_ends() */
