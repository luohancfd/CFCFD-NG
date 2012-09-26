/** \file l_misc.cxx
 * \ingroup l1d3
 * \brief Miscellaneous functions for l1d.c.
 *
 * \author PA Jacobs
 * \version ANSI Compiler version 17-Aug-95
 * \version 20-Jul-06, C++ port.
 */

/*-----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <valarray>
#include "../../../lib/util/source/useful.h"
#include "l_kernel.hh"
#include "l_tube.hh"
#include "l1d.hh"

/*=================================================================*/

/// \brief Allocate for the internal arrays of the gas slug. 
int L_alloc(struct slug_data *A)
{
    Gas_model *gmodel = get_gas_model_ptr();
    // Check for obvious errors.
    if (A == NULL) {
        printf("\nNULL data structure pointer\n");
        exit(-1);
    }
    if (A->nxdim <= 0) {
        printf("\nDeclared dimensions are too small: nxdim=%d\n", A->nxdim);
        exit(-1);
    }
    // The initial flow state needs to be kept somewhere.
    A->init_str = (struct L_flow_state *) calloc(1, sizeof(struct L_flow_state));
    // An array of cells with internal structures that need to be allocated.
    A->Cell = (struct L_cell *) calloc( A->nxdim, sizeof(struct L_cell) );
    if (A->Cell == NULL) {
        printf("\nL_alloc(): could not allocate memory.\n");
	exit( MEMORY_ERROR );
    } else {
	for ( int i = 0; i < A->nxdim; ++i ) {
	    A->Cell[i].gas = new Gas_data(gmodel);
	    A->Cell[i].ref = new Gas_data(gmodel);
	}
    }
    printf("L_alloc(): done allocating memory.\n");
    return SUCCESS;
} // end function L_alloc()

void L_free(struct slug_data *A)
{
    for ( int i = 0; i < A->nxdim; ++i ) {
	delete A->Cell[i].gas;
	delete A->Cell[i].ref;
    }
    free(A->Cell);
    delete A->init_str->gas;
    free(A->init_str);
}

/// \brief Set up min and max indices for convenience in later work.
///
/// Active cells should then be addressible as 
/// cell[ix], ixmin <= ix <= ixmax.
int L_set_index_range(struct slug_data *A)
{
    A->ixmin = A->nghost;
    A->ixmax = A->ixmin + A->nnx - 1;
    return 0;
} /* end function */


/// \brief Copy all of the gas-dynamic data from the source cell to the destination cell.
int L_copy_cell_data(struct L_cell *source,
                     struct L_cell *destination, int copy_extras)
{
    // Basic copy: just enough for the application of boundary
    // conditions to the end of the gas slugs.
    destination->gas->copy_values_from(*(source->gas));
    destination->xmid = source->xmid;
    destination->L_bar = source->L_bar;
    destination->u = source->u;
    if (copy_extras == 1) {
	// Now, for the extra items that will be needed for cell refinement.
        destination->x = source->x;
        destination->area = source->area;
        destination->pface = source->pface;
        destination->uface = source->uface;
        destination->volume = source->volume;

        destination->T_Wall = source->T_Wall;
        destination->K_over_L = source->K_over_L;
        destination->shear_stress = source->shear_stress;
        destination->heat_flux = source->heat_flux;
        destination->entropy = source->entropy;

        destination->mass = source->mass;
        destination->moment = source->moment;
        destination->Energy = source->Energy;
    }
    return 0;
} // end function L_copy_cell_data()


/// \brief Interpolate gas properties at a particular x-location.
///
/// If the x-location falls within the gas slug, the values 
/// are interpolated from the adjacent cells and the function returns 1.
/// If the location does not fall within the bounds of the 
/// gas slug, the values are set to zero and the function returns 0.
/// Look at the code to see what values correspond to what properties.
int L_interpolate_cell_data(struct slug_data *A, 
			    double xloc, struct L_cell& icell)
{
    double alpha, xmid_m1, xmid_p1, xmid, dx_plus, dx_minus;
    struct L_cell *c0, *c1;

    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();

    int found = 0;
    int ixmin = A->ixmin;
    int ixmax = A->ixmax;
    if (xloc > A->Cell[ixmin - 1].x && xloc <= A->Cell[ixmax].x) {
        // ...and then find the cell containing the x-location.
        for ( int ix = ixmin; ix <= ixmax; ++ix ) {
            if ( xloc <= A->Cell[ix].x ) {
                /* 
                 * We have found the cell 
                 * Linearly interpolate the quantities.
                 */
                xmid = 0.5 * (A->Cell[ix-1].x + A->Cell[ix].x);
                if (ix > ixmin) {
                    xmid_m1 = 0.5 * (A->Cell[ix-2].x + A->Cell[ix-1].x);
                } else {
                    xmid_m1 =
                        A->Cell[ixmin-1].x - (xmid - A->Cell[ixmin-1].x);
                }   /* end if */
                if (ix < ixmax) {
                    xmid_p1 = 0.5 * (A->Cell[ix].x + A->Cell[ix+1].x);
                } else {
                    xmid_p1 = A->Cell[ixmax].x + (A->Cell[ixmax].x - xmid);
                }   /* end if */
                dx_minus = xmid - xmid_m1;
                dx_plus = xmid_p1 - xmid;

                alpha = xloc - xmid;
                if (alpha >= 0.0) {
                    alpha /= dx_plus;
		    c0 = &( A->Cell[ix] );
                    c1 = &( A->Cell[ix+1] );
                } else {
                    alpha /= -1.0 * dx_minus;  /* fixed negative 13-Jul-02 */
		    c0 = &( A->Cell[ix] );
                    c1 = &( A->Cell[ix-1] );
                }   /* end if */
		// FIX-ME, we should delegate the blonding of the gas properties to the gas_model
		icell.gas->rho = (1.0-alpha)*c0->gas->rho + alpha * c1->gas->rho;
		icell.gas->p = (1.0-alpha)*c0->gas->p + alpha*c1->gas->p;
		icell.gas->a = (1.0-alpha)*c0->gas->a + alpha*c1->gas->a;
		for ( int isp = 0; isp < nsp; ++isp ) {
		    icell.gas->massf[isp] = (1.0-alpha) * c0->gas->massf[isp] +
			alpha*c1->gas->massf[isp];
		}
		for ( int imode = 0; imode < nmodes; ++imode ) {
		    icell.gas->e[imode] = (1.0-alpha)*c0->gas->e[imode] + alpha*c1->gas->e[imode];
		    icell.gas->T[imode] = (1.0-alpha)*c0->gas->T[imode] + alpha*c1->gas->T[imode];
		}
		icell.u = (1.0-alpha) * c0->u + alpha * c1->u;
		icell.shear_stress = (1.0-alpha)*c0->shear_stress + alpha * c1->shear_stress;
		icell.heat_flux = (1.0-alpha)*c0->heat_flux + alpha*c1->heat_flux;
		icell.entropy = (1.0-alpha)*c0->entropy + alpha*c1->entropy;
                found = 1;
                break;
            }   /* end if ... */
        }   /* end for (ix... */
    }   /* end if ... */

    return found;
} // end function L_interpolate_cell_data


/// \brief Compute the average pressure in a region near the end of the gas slug.
///
/// \param A : pointer to the gas slug data structure.
/// \param which_end : integer to indicate which end of the gas slug we want
/// \param dx : distance over which we will sample the cells
/// \returns the average of the pressures within the cells.
double L_slug_end_pressure(struct slug_data *A, int which_end, double dx)
{
    int ix, n;
    struct L_cell *c;
    double x, x0, p_avg;

#   if DEBUG >= 2
    printf( "\nGet average pressure from slug end %d over distance %g...\n",
            which_end, dx );
#   endif

    x0 = 0.0;
    p_avg = 0.0;
    /* gcc could not determine if the following code initialized x0 and p_avg */
    if ( which_end == LEFT ) {
        n = 0;
        for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
            c = &(A->Cell[ix]);
            x = c->xmid;
            if ( n == 0 ) {
                x0 = x;
                p_avg = c->gas->p;
                n = 1;
            } else {
                if ( fabs(x - x0) > dx ) break;
                p_avg += c->gas->p;
                ++n;
            }
        }   /* end for */
    } else {
        /* Assume right-hand end. */
        n = 0;
        for (ix = A->ixmax; ix >= A->ixmin; --ix) {
            c = &(A->Cell[ix]);
            x = c->xmid;
            if ( n == 0 ) {
                x0 = x;
                p_avg = c->gas->p;
                n = 1;
            } else {
                if ( fabs(x - x0) > dx ) break;
                p_avg += c->gas->p;
                ++n;
            }
        }   /* end for */
    }   /* end if */
    p_avg /= n;

#   if DEBUG >= 2
    printf( "Average pressure is: %g, over %d cells\n", p_avg, n );
#   endif
    return p_avg;
} // end function L_slug_end_pressure()

/// \brief Compute the average flow properties in a region near the end of the gas slug.
///
/// \param A : pointer to the gas slug data structure.
/// \param which_end : integer to indicate which end of the gas slug we want
/// \param dx : distance over which we will sample the cells (initially)
/// \param total_mass : total mass included in dx region
/// \param Q : store the average flow properties in here
/// \returns : success or failure
int L_slug_end_properties( struct slug_data *A,
			   int which_end,
			   double dx,
			   double * total_mass,
			   struct L_flow_state * Q)
{
    /*
     * Purpose...
     * -------
     * Compute the average pressure, temperature and velocity in a region
     * near the end of the gas slug.
     *
     * Input...
     * -----
     * *A        : pointer to the gas slug data structure.
     * which_end : integer to indicate which end of the gas slug we want
     * dx        : distance over which we will sample the cells
     * QX	 : a flow state structure for averaged properties
     *
     * Returns...
     * -------
     * average pressure, temperature and velocity
     */
    int ix, n, i;
    struct L_cell *c;
    double x, x0;

#   if DEBUG >= 2
    printf( "\nGet average pressure from slug end %d over distance %g...\n",
            which_end, dx );
#   endif

    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();

    x0 = 0.0;
    Q[0].gas->p=0.0;
    Q[0].gas->T[0]=0.0;
    Q[0].u=0.0;
    /* gcc could not determine if the following code initialized x0 and p_avg */
    if ( which_end == LEFT ) {
        n = 0;
        for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
            c = &(A->Cell[ix]);
            x = c->xmid;
            if ( n == 0 ) {
                x0 = x;
                Q[0].gas->p = c->gas->p;
		Q[0].gas->T[0] = c->gas->T[0];
		Q[0].u = c->u;
		*total_mass=c->mass;
		for (i=0; i<nsp; i++) {
		    Q[0].gas->massf[i]=c->gas->massf[i];
		}
                n = 1;
            } else {
                if ( fabs(x - x0) > dx ) break;
                Q[0].gas->p += c->gas->p;
		Q[0].gas->T[0] += c->gas->T[0];
		Q[0].u += c->u;
		*total_mass +=c->mass;
                ++n;
            }
        }   /* end for */
    } else {
        /* Assume right-hand end. */
        n = 0;
        for (ix = A->ixmax; ix >= A->ixmin; --ix) {
            c = &(A->Cell[ix]);
            x = c->xmid;
            if ( n == 0 ) {
                x0 = x;
                Q[0].gas->p = c->gas->p;
		Q[0].gas->T[0] = c->gas->T[0];
		Q[0].u = c->u;
		*total_mass=c->mass;
		for (i=0; i<nsp; i++) {
		    Q[0].gas->massf[i]=c->gas->massf[i];
		}
                n = 1;
            } else {
                if ( fabs(x - x0) > dx ) break;
                Q[0].gas->p += c->gas->p;
		Q[0].gas->T[0] += c->gas->T[0];
		Q[0].u += c->u;
		*total_mass +=c->mass;
                ++n;
            }
        }   /* end for */
    }   /* end if */
    Q[0].gas->p /= n;
    Q[0].gas->T[0] /= n;
    Q[0].u /= n;

#   if DEBUG >= 2
#   endif
    return SUCCESS;
}   /* end function L_slug_end_properties() */

/*-----------------------------------------------------------------*/

#define BLEND_PUT 0
#define BLEND_RHOUE 1

int L_blend_cells( struct L_cell *cA, struct L_cell *cB, 
		   struct L_cell *c, double alpha, int blend_type )
{
    int isp;
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();

    if ( blend_type == BLEND_PUT ) {
	c->u     = (1.0 - alpha) * cA->u     + alpha * cB->u;
	c->gas->p = (1.0 - alpha) * cA->gas->p + alpha * cB->gas->p;
	c->gas->T[0] = (1.0 - alpha) * cA->gas->T[0] + alpha * cB->gas->T[0];
	for ( isp = 0; isp <= nsp; ++isp ) {
	    c->gas->massf[isp] = (1.0 - alpha) * cA->gas->massf[isp] + 
		            alpha * cB->gas->massf[isp];
	}
	gmodel->eval_thermo_state_pT(*(c->gas));
	gmodel->eval_transport_coefficients(*(c->gas));
    } else {
	c->u       = (1.0 - alpha) * cA->u       + alpha * cB->u;
	c->gas->rho = (1.0 - alpha) * cA->gas->rho + alpha * cB->gas->rho;
	c->gas->e[0]= (1.0 - alpha) * cA->gas->e[0]+ alpha * cB->gas->e[0];
	c->gas->T[0]= (1.0 - alpha) * cA->gas->T[0]+ alpha * cB->gas->T[0];
	for ( isp = 0; isp <= nsp; ++isp ) {
	    c->gas->massf[isp] = (1.0 - alpha) * cA->gas->massf[isp] + 
		            alpha * cB->gas->massf[isp];
	}
	gmodel->eval_thermo_state_rhoe(*(c->gas));
	gmodel->eval_transport_coefficients(*(c->gas));
    }
    return 0;
} // end L_blend_cells()


int L_blend_slug_ends( struct slug_data *A, int endA, 
			  struct slug_data *B, int endB,
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
    struct L_cell *c, *cA, *cB;
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


/*-----------------------------------------------------------------*/

int L_fill_data(struct slug_data *A)
{
    /*
     * Purpose...
     * -------
     * Fill the flow field for this block with uniform
     * (initial) conditions.
     *
     * Input...
     * -----
     * *A  : pointer to the gas slug data structure.
     *
     */
    int ix;
    struct L_cell *target;
    struct L_flow_state *source;

#   if (DEBUG >= 1)
    printf("\nFill block with uniform conditions...\n");
#   endif

    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        target = &(A->Cell[ix]);
        source = A->init_str;
        target->gas->copy_values_from(*(source->gas));
        target->u = source->u;
    }   /* end for */

    return 0;
}   /* end function L_fill_data */


/*-------------------------------------------------------*/

int L_compute_areas(struct slug_data *A, TubeModel *tube)
{
    /*
     * Purpose...
     * -------
     * Compute the tube cross-section area, local loss coefficient
     * and specified wall temperature as a function of position.
     * Also, compute the cell mid-points.
     *
     * Input...
     * -----
     * *A    : pointer to gas slug data structure.
     * *tube : pointer to the tube area specification
     *
     */
    int ix, i, im1;
    double alpha;

#   if DEBUG >= 2
    printf("\nCompute areas, cell-volumes, K/L and midpoints...\n");
#   endif

    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
        /*
         * Compute appropriate interval.
         */
        im1 = (int) ((A->Cell[ix].x - tube->x1) / tube->dx);
        i = im1 + 1;
        /*
         * Now interpolate the area.
         */
        if (i <= 0) {
            A->Cell[ix].area = tube->area[0];
        } else if (i >= tube->n) {
            A->Cell[ix].area = tube->area[tube->n - 1];
        } else {
            alpha = (A->Cell[ix].x - tube->x1 - im1 * tube->dx) / tube->dx;
            A->Cell[ix].area = tube->area[im1] * (1.0 - alpha) + tube->area[i] * alpha;
        }   /* end if */
    }   /* end for */

    /*
     * Compute cell volumes and cell midpoints 
     * from the current interface areas and positions.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        A->Cell[ix].volume = 0.5 * (A->Cell[ix].area + A->Cell[ix - 1].area)
            * (A->Cell[ix].x - A->Cell[ix - 1].x);
        A->Cell[ix].xmid = 0.5 * (A->Cell[ix].x + A->Cell[ix - 1].x);
    }   /* end for */

    /*
     * Interpolate the Loss coefficients and specified wall temperatures
     * at the cell midpoints.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        /*
         * Compute appropriate interval.
         */
        im1 = (int) ((A->Cell[ix].xmid - tube->x1) / tube->dx);
        i = im1 + 1;
        /*
         * Now interpolate the loss coefficient.
         */
        if (i <= 0) {
            A->Cell[ix].K_over_L = tube->K_over_L[0];
            A->Cell[ix].T_Wall = tube->T_Wall[0];
        } else if (i >= tube->n) {
            A->Cell[ix].K_over_L = tube->K_over_L[tube->n - 1];
            A->Cell[ix].T_Wall = tube->T_Wall[tube->n - 1];
        } else {
            alpha = (A->Cell[ix].xmid - tube->x1 - im1 * tube->dx)
                / tube->dx;
            A->Cell[ix].K_over_L = tube->K_over_L[im1] * (1.0 - alpha)
                + tube->K_over_L[i] * alpha;
            A->Cell[ix].T_Wall = tube->T_Wall[im1] * (1.0 - alpha)
                + tube->T_Wall[i] * alpha;
        }   /* end if */
    }   /* end for */

    return SUCCESS;
}   /* end function L_compute_areas */

/*-----------------------------------------------------------------*/

int L_dump_cell(struct slug_data *A, int ix)
{
    /*
     * Purpose...
     * -------
     * Aid the debugging process by dumping the data for a specified
     * cell in a format that is readable.
     *
     * Input...
     * -----
     * *A    : pointer to the slug data-structure
     * ix    : cell index
     *
     */
    struct L_cell *C, *Cm1;
    int isp;

    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();

    C = &(A->Cell[ix]);
    Cm1 = &(A->Cell[ix - 1]);

    printf("-------- Cell[%d] -------:\n", ix);

    printf("Iface j-1/2: x=%e, A=%e, p=%e, u=%e\n",
           Cm1->x, Cm1->area, Cm1->pface, Cm1->uface);
    printf("Iface j+1/2: x=%e, A=%e, p=%e, u=%e\n",
           C->x, C->area, C->pface, C->uface);

    printf("Vol=%e, xmid=%e\n", C->volume, C->xmid);

    printf("C_v=%e, C_p=%e, R=%e\n", gmodel->Cv(*(C->gas)), 
	   gmodel->Cp(*(C->gas)), gmodel->R(*(C->gas)));
    for ( isp = 0; isp < nsp; ++isp ) {
        printf("f[%d]=%e ", isp, C->gas->massf[isp]);
    }
    printf("\n");
    printf("rho=%e, u=%e, e=%e\n", C->gas->rho, C->u, C->gas->e[0]);
    printf("p=%e, a=%e, T=%e, mu=%e\n", C->gas->p, C->gas->a, C->gas->T[0], C->gas->mu);
    printf("mass=%e, moment=%e, Energy=%e\n", C->mass, C->moment, C->Energy);
    printf("Q_m=%e, Q_mom=%e, Q_E=%e\n", C->Q_m, C->Q_mom, C->Q_E);

    printf("--------------------------------------\n");

    return 0;
}   /* end of L_dump_data() */

/*---------------------------------------------------------------*/

int maximum_p(struct slug_data *A, double *p_max, double *x_max)
{
    /* Purpose...
     * -------
     * For the given slug of gas, find the maximum pressure and
     * its location.
     * 06-Oct-92  -- smooth pressures before search
     * 
     * Input...
     * -----
     *  *A   : pointer to the slug_data structure
     *
     * Output...
     * ------
     *  *p_max  : pointer to the value of maximum pressure
     *  *x_max  : pointer to the x-location of the max pressure
     *
     */
    double xx, pp, p;
    int ix, something_done;

    pp = 0.0;
    xx = 0.0;
    something_done = 0;

    for (ix = A->ixmin + 1; ix <= A->ixmax - 1; ++ix) {
        p = (A->Cell[ix - 1].gas->p + A->Cell[ix].gas->p + 
	     A->Cell[ix + 1].gas->p) / 3.0;
        if (p > pp) {
            pp = p;
            xx = A->Cell[ix].xmid;
        }   /* end if */
        something_done = 1;
    }   /* end for */

    if ( !something_done ) {
        /* Arbitrary choice as there are not enough cells, anyway. */
        pp = A->Cell[A->ixmin].gas->p;
        xx = A->Cell[A->ixmin].xmid;
    }

    *p_max = pp;
    *x_max = xx;

    return (0);
}   /* end function maximum_p */

/*--------------------------------------------------------------*/

int total_energy(struct slug_data *A, double *E_tot)
{
    /* Purpose...
     * -------
     * For the given slug of gas, compute the total energy
     * Added 20-Apr-93
     * 
     * Input...
     * -----
     *  *A   : pointer to the slug_data structure
     *
     * Output...
     * ------
     *  *E_tot  : pointer to the value of total energy
     *
     */
    double sum;
    int ix;

    sum = 0.0;

    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        sum += A->Cell[ix].mass * 
               (A->Cell[ix].gas->e[0] + 0.5 * A->Cell[ix].u * A->Cell[ix].u);
    }   /* end for */

    *E_tot = sum;

    return 0;
}   /* end function total_energy */

/*------------------------------------------------------------*/

double L_get_dt_plot(SimulationData *SD) 
{
    int i;
    double dt_plot;
    dt_plot = SD->dt_plot[0];
    for ( i = 0; i < SD->n_dt_plot; ++i ) {
        if ( SD->sim_time > SD->t_change[i] ) {
            dt_plot = SD->dt_plot[i];
        }
    }
    return dt_plot;
}   /* end L_get_dt_plot() */

double L_get_dt_history(SimulationData *SD) 
{
    int i;
    double dt_his;
    dt_his = SD->dt_his[0];
    for ( i = 0; i < SD->n_dt_plot; ++i ) {
        if ( SD->sim_time > SD->t_change[i] ) {
            dt_his = SD->dt_his[i];
        }
    }
    return dt_his;
}   /* end L_get_dt_history() */

/*============== end of l_misc.c ================*/
