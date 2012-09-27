/** \file l_adapt.cxx
 * \ingroup l1d2
 * \brief Cell-adaption routines for l1d.c.
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../../../lib/util/source/useful.h"
#include "l1d.hh"
#include "l_kernel.hh"
#include "l_cell.hh"
#include "l_adapt.hh"

/*-----------------------------------------------------------------*/

#define LEAVE_AS_IS 0
#define FUSE_CELL   1
#define SPLIT_CELL  2

/*-----------------------------------------------------------------*/

int test_for_expansion(double dudx[], int indicator[], int ixmin, int ixmax)
{
    /* Detect shocks of sufficient strength.
     * Weak jumps may be caused by other disturbances such as
     * cell merging or splitting.
     */
    int ix;
    double dudx_max;

    /* Locate the largest expansion rate. */
    dudx_max = 1.0;
    for (ix = ixmin; ix <= ixmax; ++ix) {
        if (dudx[ix] > dudx_max)
            dudx_max = dudx[ix];
    }   /* end for */

    for (ix = ixmin; ix <= ixmax; ++ix) {
        if (dudx[ix] > 0.1 * dudx_max) {
            indicator[ix] = 1;
        } else {
            indicator[ix] = 0;
        }   /* end if */
    }   /* end for */

    return 0;
}   /* end function test_for_expansion() */

/*-----------------------------------------------------------------*/

int test_for_shock(double dudx[], int indicator[], int ixmin, int ixmax)
{
    /* Detect shocks of sufficient strength.
     * Weak jumps may be caused by other disturbances such as
     * cell merging or splitting.
     */
    int ix, result;
    double dudxL, dudx_i, dudxR;

    for (ix = ixmin; ix <= ixmax; ++ix) {
        if (ix == ixmin) {
            dudxL = dudx[ix];
            dudx_i = dudx[ix];
            dudxR = dudx[ix + 1];
        } else if (ix == ixmax) {
            dudxL = dudx[ix - 1];
            dudx_i = dudx[ix];
            dudxR = dudx[ix + 1];
        } else {
            dudxL = dudx[ix - 1];
            dudx_i = dudx[ix];
            dudxR = dudx[ix];
        }   /* end if */

        result = 0;
        if (dudx_i < 0.0) {
            /* We are in a compression region; is it strong? */
            if (dudxL < 1.1 * dudx_i || dudxL > 0.9 * dudx_i)
                result = 1;
            if (dudxR < 1.1 * dudx_i || dudxR > 0.9 * dudx_i)
                result = 1;
        }
        /* end if */
        indicator[ix] = result;
    }   /* end for */

    return 0;
}   /* end function test_for_shock() */

/*-----------------------------------------------------------------*/

int shock_is_near(int shk[], int ix, int ixmin, int ixmax)
{
    /* Look for shocks nearby...
     */
    int sum;

    if (ix == ixmin) {
        sum = shk[ix] + shk[ix + 1] + shk[ix + 2];
    } else if (ix == ixmin + 1) {
        sum = shk[ix - 1] + shk[ix] + shk[ix + 1] + shk[ix + 2];
    } else if (ix == ixmax - 1) {
        sum = shk[ix - 2] + shk[ix - 1] + shk[ix] + shk[ix + 1];
    } else if (ix == ixmax) {
        sum = shk[ix - 2] + shk[ix - 1] + shk[ix];
    } else {
        sum = shk[ix - 2] + shk[ix - 1] + shk[ix] + shk[ix + 1] + shk[ix + 2];
    }   /* end if */

    return (sum >= 1);
}   /* end function shock_is_near */

/*-----------------------------------------------------------------*/

int near_but_not_in_shock(int shk[], int ix, int ixmin, int ixmax)
{
    /* 
     * Indicate if a shock is within a couple of cells but
     * that this cell is not in the shock.
     */
    int result;
    result = !shk[ix] && shock_is_near(shk, ix, ixmax, ixmin);
    return result;
}   /* end function near_but_not_in_shock */

/*-----------------------------------------------------------------*/

int error_indicator(double rho[], int indicator[], int ixmin, int ixmax)
{
    /* Detect cells with significant error in density.
     * This is the indicator used by Sun & Takayama
     * J. Comput. Phys Vol. 150 pp143-180. (1999)
     */
    const double alpha = 0.03;
    const double eps_refine = 0.6;
    const double eps_coarsen = 0.05;

    double phi_i, rho_im1, rho_i, rho_ip1;
    int ix, result;

    for (ix = ixmin; ix <= ixmax; ++ix) {
        if (ix == ixmin) {
            rho_im1 = rho[ix];
            rho_i = rho[ix];
            rho_ip1 = rho[ix + 1];
        } else if (ix == ixmax) {
            rho_im1 = rho[ix - 1];
            rho_i = rho[ix];
            rho_ip1 = rho[ix + 1];
        } else {
            rho_im1 = rho[ix - 1];
            rho_i = rho[ix];
            rho_ip1 = rho[ix];
        }   /* end if */

        phi_i = fabs(rho_im1 - 2.0 * rho_i + rho_ip1) /
            alpha * rho_i + fabs(rho_ip1 - rho_im1);

        result = 0;
        if (phi_i > eps_refine) {
            result = SPLIT_CELL;
        } else if (phi_i < eps_coarsen) {
            result = FUSE_CELL;
        } else {
            result = LEAVE_AS_IS;
        }   /* end if */
        indicator[ix] = result;
    }   /* end for */

    return 0;
}   /* end function error_indicator() */

/*-----------------------------------------------------------------*/

int L_adapt_cells(GasSlug* A)
{
    /*
     * Returns 0 if OK; 1 otherwise.
     */
    int ix, ia, ib, B_nnx;
    static std::vector<LCell> B_Cell;
    LCell *ci, *cim1, *cip1;
    static int too_large[NDIM], too_small[NDIM];
    static int exp_indicator[NDIM], shock_indicator[NDIM];
    static int density_indicator[NDIM], decision[NDIM];
    int should_be_refined, should_be_coarsened;
    static double dudx[NDIM], rho[NDIM];
    double dx, dx_min, dx_max, dx_refine, dx_coarsen;
    Gas_model *gmodel = get_gas_model_ptr();

    if ( B_Cell.size() == 0 ) {
	// If one cell is not filled out, assume the rest are not.
	for ( ix = 0; ix < NDIM; ++ix ) {
	    B_Cell.push_back(LCell(gmodel));
	}
	printf( "L_adapt_cells(): first-time: work arrays filled out.\n" );
    }

#   define SMOOTH_SPLIT 1
    /* Set the above parameter to 0 to get the old, rough splitting. */

    if (A->adaptive == ADAPT_NONE)
        return 0;

    dx_min = A->dxmin;
    dx_max = A->dxmax;
    dx_refine = dx_min * 4.0;
    dx_coarsen = dx_max / 4.0;

    for (ix = 0; ix < NDIM; ++ix) {
        too_small[ix] = 0;
        too_large[ix] = 0;
        decision[ix] = 0;
        dudx[ix] = 0.0;
    }   /* end for */

    /*
     * Identify large & small cells.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        dx = A->Cell[ix].x - A->Cell[ix - 1].x;
        if (fabs(dx) < dx_min) {
            too_small[ix] = 1;
        }
        if (fabs(dx) > dx_max) {
            too_large[ix] = 1;
        }
    }   /* end for */

    /*
     * Identify cells in or near shocks and expansion fans.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        dx = A->Cell[ix].x - A->Cell[ix - 1].x;
        dudx[ix] = (A->Cell[ix].u - A->Cell[ix - 1].u) / dx;
        rho[ix] = A->Cell[ix].gas->rho;
    }   /* end for */

    test_for_shock(dudx, shock_indicator, A->ixmin, A->ixmax);
    test_for_expansion(dudx, exp_indicator, A->ixmin, A->ixmax);
    error_indicator(rho, density_indicator, A->ixmin, A->ixmax);

    /*
     * Decide what to do with each cell.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        dx = A->Cell[ix].x - A->Cell[ix - 1].x;

        should_be_refined = (too_large[ix] ||
            (density_indicator[ix] == SPLIT_CELL
              && fabs(dx) > dx_refine && exp_indicator[ix])
            );

        should_be_coarsened = ( (too_small[ix] ||
            (density_indicator[ix] == FUSE_CELL && fabs(dx) < dx_coarsen))
            && !shock_is_near(shock_indicator, ix, A->ixmin, A->ixmax)
            );

        decision[ix] = LEAVE_AS_IS;
        if (should_be_refined && (A->adaptive != ADAPT_FUSE_ONLY)) {
            decision[ix] = SPLIT_CELL;
        } else if (should_be_coarsened && (A->adaptive != ADAPT_SPLIT_ONLY)) {
            decision[ix] = FUSE_CELL;
        }   /* end if */
    }   /* end for */

    /*
     * Copy, with merge or split, to the temporary array.
     */
    ia = A->ixmin;
    ib = A->ixmin;
    while (ia < A->ixmax) {
	ci = &(A->Cell[ia]);
	cim1 = &(A->Cell[ia-1]);
	cip1 = &(A->Cell[ia+1]);
        if (decision[ia] == SPLIT_CELL) {
#           if SMOOTH_SPLIT == 1
	    L_split_cell_data_smoothly( ci, cim1, cip1, &(B_Cell[ib]), 
					&(B_Cell[ib + 1]) );
#           else
            L_split_cell_data_roughly( ci, cim1->x, &(B_Cell[ib]), 
				       &(B_Cell[ib + 1]) );
#           endif
            ia += 1;
            ib += 2;
        } else if (decision[ia] == FUSE_CELL || decision[ia + 1] == FUSE_CELL) {
            L_fuse_cell_data( ci, cip1, &(B_Cell[ib]) );
            ia += 2;
            ib += 1;
        } else {
            L_copy_cell_data( ci, &(B_Cell[ib]), 1 );
            ia += 1;
            ib += 1;
        }   /* end if */
    }   /* end while */

    /* If we have one left over, copy it or split it. */
    if (ia == A->ixmax) {
	ci = &(A->Cell[ia]);
	cim1 = &(A->Cell[ia-1]);
	cip1 = &(A->Cell[ia+1]);
        if (decision[ia] == SPLIT_CELL) {
#           if SMOOTH_SPLIT == 1
	    L_split_cell_data_smoothly( ci, cim1, cip1, &(B_Cell[ib]), 
					&(B_Cell[ib + 1]) );
#           else
            L_split_cell_data_roughly( ci, cim1->x, &(B_Cell[ib]), 
				       &(B_Cell[ib + 1]) );
#           endif
            ia += 1;
            ib += 2;
        } else {
            L_copy_cell_data( ci, &(B_Cell[ib]), 1 );
            ia += 1;
            ib += 1;
        }   /* end if */
    }   /* end if */

    /*
     * Copy back into the slug data array.
     */
    B_nnx = ib - A->ixmin;
    if (B_nnx <= A->nxmax) {
        A->nnx = B_nnx;
	A->set_index_range();
        for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
            L_copy_cell_data(&(B_Cell[ix]), &(A->Cell[ix]), 1);
        }   /* end for */
        if ( A->decode_conserved() != 0 ) {
	    printf( "L_adapt_cells(): Failure decoding conserved quantities.\n" );
	    return 1;
	}
    } else {
        printf("L_adapt_cells(): Could not adapt cells, not enough space.\n");
        printf("A.nxmax = %d, B_nnx = %d\n", A->nxmax, B_nnx);
    }   /* end if */

    return 0;
}   /* end function L_adapt_cells() */

/*-----------------------------------------------------------------*/

int L_fuse_cell_data(LCell *source1, LCell *source2, LCell *destination)
{
    /*
     * Fuse the gas-dynamic data from the two source cells
     * into the one destination cell.
     * It is assumed that source1 is to the left of and
     * adjacent to source2.
     */

    destination->xmid = 0.5 * (source1->xmid + source2->xmid);
    destination->L_bar = 0.5 * (source1->L_bar + source2->L_bar);

    /* 
     * This averaging probably doesn't matter so much.
     * It will be written over by during a later call to L_decode_conserved().
     * But, just in case something is call upon in the mean time...
     */
    destination->gas->copy_values_from(*(source2->gas));
    destination->gas->rho = 0.5 * (source1->gas->rho + source2->gas->rho);
    destination->u       = 0.5 * (source1->u + source2->u);
    destination->gas->e[0]= 0.5 * (source1->gas->e[0] + source2->gas->e[0]);
    destination->gas->p   = 0.5 * (source1->gas->p + source2->gas->p);
    destination->gas->a   = 0.5 * (source1->gas->a + source2->gas->a);
    destination->gas->T[0]= 0.5 * (source1->gas->T[0] + source2->gas->T[0]);

    /*
     * The unified cell is bounded on the right by source2's 
     * interface; copy the interface values from there.
     * The interface between the two is forgotten.
     */
    destination->x = source2->x;
    destination->area = source2->area;
    destination->pface = source2->pface;
    destination->uface = source2->uface;
    destination->volume = source1->volume + source2->volume;

    /*
     * Combine the conserved quantities for both original cells.
     * This is the real merge.
     */
    destination->mass = source1->mass + source2->mass;
    destination->moment = source1->moment + source2->moment;
    destination->Energy = source1->Energy + source2->Energy;

    return 0;
}   /* end function L_fuse_cell_data() */

/*-----------------------------------------------------------------*/

int L_split_cell_data_roughly(LCell *source, double xL, LCell *dest1, LCell *dest2)
{
    /*
     * Split the gas-dynamic data from the one source cell
     * into the two destination cells.
     * It is assumed that dest1 is to the left of and
     * adjacent to dest2.
     */

    dest1->x = 0.5 * (xL + source->x);
    dest2->x = source->x;
    dest1->xmid = 0.5 * (xL + dest1->x);
    dest2->xmid = 0.5 * (dest1->x + dest2->x);

    dest1->area = source->area;
    dest1->pface = source->pface;
    dest1->uface = source->uface;
    dest1->volume = 0.5 * source->volume;

    dest2->area = source->area;
    dest2->pface = source->pface;
    dest2->uface = source->uface;
    dest2->volume = 0.5 * source->volume;

    dest1->L_bar = source->L_bar;
    dest2->L_bar = source->L_bar;

    /*
     * Copy the same gas properties into both cells.
     */
    dest1->gas->copy_values_from(*(source->gas));
    dest2->gas->copy_values_from(*(source->gas));
    /*
     * Copy the same intensive properties into both cells.
     */
    dest1->u = source->u;
    dest2->u = source->u;

    /*
     * Split the extensive properties evenly between the two cells.
     * This may be OK for regions with mild gradients but may be
     * incorrect for regions with large flow gradients or with large
     * area gradients.
     */
    dest1->mass = 0.5 * source->mass;
    dest1->moment = 0.5 * source->moment;
    dest1->Energy = 0.5 * source->Energy;

    dest2->mass = 0.5 * source->mass;
    dest2->moment = 0.5 * source->moment;
    dest2->Energy = 0.5 * source->Energy;

    return 0;
}   /* end function L_split_cell_data_roughly() */

/*-----------------------------------------------------------------*/

int L_split_cell_data_smoothly(LCell *ci,   /* the cell to be split */
			       LCell *cim1, /* the cell to left     */
			       LCell *cip1, /* the cell to right    */
			       LCell *ca,   /* left new cell        */ 
			       LCell *cb )  /* right new cell       */
{
    /*
     * Split the gas-dynamic data from the one source cell
     * into the two destination cells.
     * It is assumed that ca is to the left of and adjacent to cb.
     * See workbook 2001/1 for details.
     */
    double msplit, esplit;   /* splitting parameters for mass and energy */
    double xs, As;           /* new interface position and area          */
    double s1, s2, s_rho, s_e, dxi;
    UNUSED_VARIABLE(s_e);

#   define DEBUG_NEW_SPLIT 0
#   if DEBUG_NEW_SPLIT == 1
    printf( "Begin L_split_cell_data_smoothly()\n" );
#   endif

    /*
     * Geometry
     */
    dxi = ci->x - cim1->x;
    xs = 0.5 * (ci->x + cim1->x);
    ca->x = xs;
    cb->x = ci->x;
    ca->xmid = 0.5 * (cim1->x + ca->x);
    cb->xmid = 0.5 * (ca->x + cb->x);

    As = 0.5 * (ci->area + cim1->area);  /* assume linear variation */
    ca->area = As;
    cb->area = ci->area;

    ca->volume = 0.5 * (cim1->area + As) * (xs - cim1->x);
    cb->volume = 0.5 * (As + ci->area) * (ci->x - xs);
    if ( (ca->volume + cb->volume - ci->volume) / ci->volume > 0.01 ) {
	printf( "L_split_cell_data_smoothly(): areas don't add up\n" );
	printf( "Va=%g, Vb=%g, V=%g\n", ca->volume, cb->volume, ci->volume );
    }

    /*
     * We had better copy the following values, 
     * in case they are later used.
     */
    ca->pface = ci->pface;
    ca->uface = ci->uface;
    cb->pface = ci->pface;
    cb->uface = ci->uface;

    ca->L_bar = ci->L_bar;
    cb->L_bar = ci->L_bar;

    /*
     * Copy the same gas properties into both cells.
     * This will make sure that we have the mass fractions copied.
     */
    ca->gas->copy_values_from(*(ci->gas));
    cb->gas->copy_values_from(*(ci->gas));
    /*
     * Copy the same intensive properties into both cells.
     */
    ca->u = ci->u;
    cb->u = ci->u;

    /*
     * Split the mass between the two new cells such that
     * the local slope in density is maintained.
     */
    s1 = (ci->gas->rho - cim1->gas->rho) / (ci->x - cim1->x);
    s2 = (cip1->gas->rho - ci->gas->rho) / (cip1->x - ci->x);
#   if 0
    /* A geometric mean of the slopes will be conservative. */
    if ( s1 * s2 > 0.0 ) 
	s_rho = 2.0 * s1 * s2 / (s1 + s2);
    else
	s_rho = 0.0;
#   else
    /* An arithmetic mean might be more responsive. */
    s_rho = 0.5 * (s1 + s2);
#   endif
    msplit = (s_rho * dxi * 0.5 / ci->mass + 1.0 / ca->volume ) /
	( 1.0 / ca->volume + 1.0 / cb->volume );
    ca->mass = (1.0 - msplit) * ci->mass;
    cb->mass = msplit * ci->mass;

#   if DEBUG_NEW_SPLIT == 1
    printf( "s1=%g, s2=%g, s_rho=%g, msplit=%g\n", 
	    s1, s2, s_rho, msplit );
#   endif

    /*
     * The momentum can now be reconstructed.
     */
    ca->moment = ca->mass * ca->u;
    cb->moment = cb->mass * cb->u;

    /*
     * For the moment, have equal internal energy.  --- FIX ME ---
     * Total energy is then just split as for mass.
     */
    ca->gas->e[0] = ci->gas->e[0];
    cb->gas->e[0] = ci->gas->e[0];
    esplit = msplit;
    ca->Energy = (1.0 - esplit) * ci->Energy;
    cb->Energy = esplit * ci->Energy;

#   if DEBUG_NEW_SPLIT == 1
    printf( "End L_split_cell_data_smoothly()\n" );
#   endif
    return 0;
}   /* end function L_split_cell_data_smoothly() */

/*-----------------------------------------------------------------*/
