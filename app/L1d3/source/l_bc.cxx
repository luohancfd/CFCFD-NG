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
#include "l_misc.hh"

/*=================================================================*/

int L_bc_left_velocity(struct slug_data *A, double v)
{
    /* 
     * Apply a specified velocity to the left boundary.
     */
    int ix;
    struct L_cell *src, *dest;
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

/*------------------------------------------------------------*/

int L_bc_left_reflect(struct slug_data *A)
{
    /*
     * Apply a reflective bc to the left boundary.
     */
    int ix;
    struct L_cell *src, *dest;
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

/*------------------------------------------------------------*/

int L_bc_left_free(struct slug_data *A)
{
    /*
     * Apply a free/supersonic outflow bc to the left boundary.
     */
    int ix;
    struct L_cell *src, *dest;
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

/*------------------------------------------------------------*/

int L_bc_right_velocity(struct slug_data *A, double v)
{
    /*
     * Apply a specified velocity to the right boundary.
     */
    int ix;
    struct L_cell *src, *dest;
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


/*------------------------------------------------------------*/

int L_bc_right_reflect(struct slug_data *A)
{
    /*
     * Apply a reflective bc to the right boundary.
     */
    int ix;
    struct L_cell *src, *dest;
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

/*------------------------------------------------------------*/

int L_bc_right_free(struct slug_data *A)
{
    /* 
     * Apply a free/supersonic-outflow bc to the right boundary.
     */
    int ix;
    struct L_cell *src, *dest;
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

/*------------------------------------------------------------*/

int L_exchange_bc_data(struct slug_data *A, struct slug_data *B)
{
    /*
     * Copy the end data from the left slug (A) to the  right slug (B).
     *       ia-1 ia    |  ia+1 ia+2
     *        V    V    |   ^    ^
     *       ib-2 ib-1  |  ib   ib+1
     */
    int ia, ib;
    struct L_cell *src, *dest;

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

/*=================================================================*/
