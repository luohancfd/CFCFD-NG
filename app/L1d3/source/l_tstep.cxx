/** \file l_tstep.cxx
 * \ingroup l1d3
 * \brief Time-stepping routines for l1d.cxx.
 *
 * \author PA Jacobs
 */

/*-----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"
#include "l_kernel.hh"
#include "l1d.hh"
#include "l_tstep.hh"
#include "l_rivp.hh"
#include "l_misc.hh"
#include "../../../lib/nm/source/qd_power.h"
#include "../../../lib/nm/source/qd_log10.h"

/*-----------------------------------------------------------------*/

/** \brief mim-mod limiter */
#define MIN_MOD(c,d)   ( ((c)*(d) <= 0.0) ? 0.0 : SMALLEST(c,d) )

/*=================================================================*/


/*
 * -------------------------------
 * Conserved variable manipulation
 * -------------------------------
 */

int L_encode_conserved_for_one_cell(struct L_cell *c)
{
    /*
     * Purpose...
     * -------
     * Compute the conserved quantities in the specified cell
     * from the primary variables.
     */

    /* Mass.  Assuming that the cell volume is correct. */
    c->mass = c->gas->rho * c->volume;

    /* X-momentum. */
    c->moment = c->mass * c->u;

    /*
     * Total Energy = mass *
     *     (specific internal energy + kinetic energy/unit mass).
     */
    c->Energy = c->mass * (c->gas->e[0] + 0.5 * (c->u * c->u));

    return 0;
}   /* end function L_encode_conserved_for_one_cell() */


int L_encode_conserved(struct slug_data *A)
{
    /*
     * Purpose...
     * -------
     * Compute the conserved quantities in each cell from the
     * primary variables.
     */
    int ix;

#   if (DEBUG >= 1)
    printf("\nEncode conserved quantities...\n");
#   endif

    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
	L_encode_conserved_for_one_cell( &(A->Cell[ix]) );
    }

    return 0;
}   /* end function L_encode_conserved */

/*-----------------------------------------------------------------*/

int L_decode_conserved(struct slug_data *A)
{
    /*
     * Compute the primary variables from the conserved quantities
     * in each cell.
     *
     * Input...
     * *A         : a pointer to the slug_data structure
     * Output...
     * Returns 0 if all is OK, 1 otherwise.
     */
    int    ix;
    double ke, gam;
    struct L_cell *c;
    Gas_model *gmodel = get_gas_model_ptr();

#   if DEBUG >= 2
    printf("Decode conserved quantities...\n");
#   endif

    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        c = &( A->Cell[ix] );
        c->gas->rho = c->mass / c->volume;
        c->u = c->moment / c->mass;
        ke = 0.5 * (c->u * c->u);
        c->gas->e[0] = c->Energy / c->mass - ke;
        if ( c->gas->rho <= 0.0 || c->mass <= 0.0 || fabs(c->u) > 1.0e5 ) {
            printf( "L_decode_conserved: Bad value for density, mass or velocity\n");
	    printf( "    rho=%g, mass=%g, u=%g in cell %d\n", c->gas->rho, c->mass, c->u, ix );
	    printf( "    x.left=%g, x.right=%g\n", A->Cell[ix-1].x, c->x );
	    printf( "    mass=%g, velocity=%g, volume=%g\n", c->mass, c->u, c->volume );
	    printf( "    momemtum=%g, total Energy=%g e[0]=%g\n", c->moment, c->Energy, c->gas->e[0] );
	    c->gas->print_values();
            exit(BAD_CELLS_ERROR);
	}

        /* Fill out the other variables; P, T, a. */
	gmodel->eval_thermo_state_rhoe(*(c->gas));
	gmodel->eval_transport_coefficients(*(c->gas));

        /* Entropy referenced to 1 atm and 300K */
        gam = gmodel->Cp(*(c->gas)) / gmodel->Cv(*(c->gas));
        c->entropy = gmodel->Cv(*(c->gas)) * 
            log(pow(c->gas->T[0] / 300.0, gam) *
                pow(c->gas->p / 101.3e3, (1.0 - gam)));
    }   /* ix loop */

    return 0;
}   /* end function L_decode_conserved */


/*-----------------------------------------------------------------*/

int L_set_chemistry_timestep(struct slug_data *A, double dt)
{
    struct L_cell *c;
    for ( int ix = A->ixmin; ix <= A->ixmax; ++ix ) {
        c = &( A->Cell[ix] );
	c->dt_chem = dt;
    }
    return 0;
}

int L_chemical_increment(struct slug_data *A, double dt)
{
    /* Use Rowan's chemistry module to update the species mass fractions. */
    int    ix, flag;
    struct L_cell *c;
    Gas_model *gmodel = get_gas_model_ptr();
    Reaction_update *rupdate = get_reaction_update_ptr();

    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        c = &( A->Cell[ix] );
	flag = rupdate->update_state(*(c->gas), dt, c->dt_chem, gmodel);
	if ( flag != 0 ) {
	    cout << "chemical update failed" << endl
		 << "    for cell " << ix << " at x=" << c->xmid << endl;
	}
	// The update only changes mass fractions, we need to impose
	// a thermodynamic constraint based on a call to the equation
	// of state.
	gmodel->eval_thermo_state_rhoe(*(c->gas));
	// Ensure viscous properties are up-to-date.
	gmodel->eval_transport_coefficients(*(c->gas));
    } /* end ix loop */

    return 0;
} /* end function L_chemical_increment */

/*----------------------------------------------------------------*/

/*
 * -----------------
 * Production Terms.
 * -----------------
 */

double f_darcy_weisbach( double Re, double D )
{
    /*
     * Darcy-Weisbach friction factor for fully-developed pipe flow.
     * Re is based on pipe diameter.
     */
    double f, lgtmp;

    /* size of wall roughness elements in metres */
    double eps = 0.025e-3;

    if (Re < 10.0) {
	/* A reasonable limit for very low speeds. */
	f = 6.4;
    } else if (Re < 2000.0) {
	/* Laminar regime. */
	f = 64.0 / Re;
    } else if (Re < 4000.0) {
	/* Transition regime. */
	f = 0.032 / pow(Re/2000.0, 0.3187);
    } else {
	/* Fully turbulent */
#       define  SMOOTH  1
#       if (SMOOTH == 1)
	lgtmp = 1.327359 - 0.9 * log10(Re);
	UNUSED_VARIABLE(D);
	UNUSED_VARIABLE(eps);
#       else
	lgtmp = log10(21.25 * pow(Re, -0.9) + eps / D);
#       endif
	f = 1.14 - 2.0 * lgtmp;
	f = 1.0 / (f * f);
    }

    return f;
}

double f_flat_plate( double Re )
{
    /* 
     * Flat-plate friction factor with Re based on length along plate.
     * See Holman (1986) Heat Transfer, Eq 5.125 to 5.127
     */
    double f;

    if (Re < 1.0e4) { 
	/* somewhat arbitrary lower Re limit to stop 
	 * possibility of near-infinite friction and heat flux */
	f = 8 * 0.332 / sqrt(10000.0);
    } else if (Re < 5.0e5) { /* Laminar regime */
	f = 8 * 0.332 / sqrt(Re); 
    } else if (Re < 1.0e7) { /* lower Re turbulent regime. */
	f = 8 * 0.0296 * pow(Re, -0.2);
    } else { /* high Re tubulent regime. */
	f = 8 * 0.185 * pow(log10(Re), -2.584);
    }

    return f;
}


int L_source_vector(struct slug_data *A)
{
    /*
     * Purpose...
     * -------
     * Compute the components of the source vector, Q.
     * This vector is used to include viscous losses and
     * heat transfer and mass sources/sinks.
     * 
     * Version...
     * -------
     * 06-Oct-93 Added a radiation term
     * 14-Nov-94 Included the BL mass-loss terms
     * 12-Feb-03 Major clean-up to include David Buttsworth's work.
     *
     */
    int ix;
    double lambda, f, Re_D, Re_L, D, abs_u, area, M;
    double T_wall, T_aw, F_LOSS, T_ref, St;
    double omega, tau0, h, length;
    double mypi, w_dot, q;
    double sigma, alpha, Tgas, q_rad;
    struct L_cell *cell;    /* pointer to the local cell */
    double r2r1, m_loss, mom_loss, E_loss, kt, kl;
    double gam;
    Gas_model *gmodel = get_gas_model_ptr();

#   if DEBUG >= 3
    printf("\nL_source_vector(): Begin...\n");
#   endif

    mypi = 3.1415927;

    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        cell = &(A->Cell[ix]);
        /*
         * Initialize the source terms to zero.
         * Each physical effect modelled will add or subtract 
	 * something from these quantities.
         */
        cell->Q_m = 0.0;
        cell->Q_mom = 0.0;
        cell->Q_E = 0.0;

	/*
	 * The following variables are used to record the magnitude
	 * of the viscous effects -- clear them in case they don't get
	 * set on this pass..
	 */
	cell->shear_stress = 0.0;
	cell->heat_flux = 0.0;

        if ( A->viscous_effects == INVISCID ) {
	    continue;  /* Nothing further to be done for this cell. */
	}

	if ( L_get_case_id() == DRUMMOND_WITH_M4_NOZZLE ) {
	    /*
	     * Treat nozzle flow as inviscid because 
	     * turbulent pipe flow model or other shocked gas models
	     * are unrealistic for losses in the (assumed steady) nozzle flow.
	     */
	    if (cell->x >= 0.049) { /* start of nozzle contraction */
		continue;  /* Nothing further to be done for this cell. */
	    }
	}

	/*
	 * Prepare to compute the viscous effects / losses.
	 */
        gam = gmodel->Cp(*(cell->gas)) / gmodel->Cv(*(cell->gas));

	/* Local tube geometry */
	area = 0.5 * (cell->area + A->Cell[ix - 1].area);
	D = 2.0 * sqrt(area / mypi); /* effective diameter */
	length = cell->x - A->Cell[ix - 1].x; /* between interfaces */
#       if DEBUG >= 3
	printf("L_source_vector(): cell[%d] -------\n", ix);
	printf("                   area=%e, D=%e, L=%e\n", area, D, length);
	printf("                   x(j-1/2)=%e, x(j+1/2)=%e\n",
	       A->Cell[ix - 1].x, cell->x);
#       endif

	/* Local gas speed */
	abs_u = fabs(cell->u);
	/* Recovery factor -- assume turbulent */
	double Prandtl = 0.75; // a constant value for the moment -- FIX-ME
	omega = pow(Prandtl, 0.333);
	/* Local Mach number */
	M = abs_u / cell->gas->a;

#       define APPLY_COMPRESSIBILITY_FACTOR 1
	if ( APPLY_COMPRESSIBILITY_FACTOR == 1 ) {
	    lambda = 1.0 + (gam - 1.0) * 0.5 * omega * M * M;
	} else {
	    lambda = 1.0;
	}

#       if DEBUG >= 3
	printf("L_source_vector(): abs_u=%e, omega=%e\n", abs_u, omega);
	printf("                   M=%e, lambda=%e\n", M, lambda);
#       endif

	/* Adiabatic wall temperature */
	T_aw = lambda * cell->gas->T[0];

	/* Local wall temperature */
	if (A->adiabatic == 1) {
	    T_wall = T_aw;
	} else if (A->adiabatic == 2) {
	    /*
	     * John Hunter's suggestion.  Heat transfer will tend to 
	     * increase the wall temperature to something close to the
	     * core temperature of the gas.
	     */
	    T_wall = cell->gas->T[0] - 400.0;
	    if (T_wall < cell->T_Wall)
		T_wall = cell->T_Wall;
	} else {
	    T_wall = cell->T_Wall;
	}

	/* 
	 * Transport properties based on Eckert reference conditions.
	 */
	T_ref = cell->gas->T[0] + 0.5 * (T_wall - cell->gas->T[0]) +
	    0.22 * (T_aw - cell->gas->T[0]);
	cell->ref->copy_values_from(*(cell->gas));
	cell->ref->T[0] = T_ref;
	cell->ref->rho = cell->gas->rho * cell->gas->T[0] / T_ref;
	cell->ref->p = cell->ref->rho * gmodel->R(*(cell->ref)) * cell->ref->T[0];
	gmodel->eval_thermo_state_pT(*(cell->ref));
	gmodel->eval_transport_coefficients(*(cell->ref));
#       if DEBUG >= 3
	printf( "L_source_vector(): T=%e, T_aw=%e, T_wall=%e, T_ref=%e, mu=%e\n",
		cell->gas->T[0], T_aw, T_wall, cell->ref->T[0], cell->ref->mu);
#       endif

	/* Local Reynolds number based on diameter and reference conditions. */
	Re_D = cell->ref->rho * D * abs_u / cell->ref->mu;
	/* use distance moved as reference distance for Re */
	Re_L = cell->ref->rho * cell->L_bar * abs_u / cell->ref->mu;
#       if DEBUG >= 3
	printf( "L_source_vector(): M=%e, lambda=%e, Re_D=%e, Re_L=%e\n",
		M, lambda, Re_D, Re_L );
#       endif

	if ( ( A->viscous_effects == VISCOUS_LOSS_FLAT_PLATE_F || 
	       A->viscous_effects == DRB_MASS_LOSS_FLAT_PLATE_F )
	     && cell->L_bar > 0.0 ) {
	    /*
	     * Friction and heat flux based on a flat plate calculation 
	     * No compressibility correction apart from
	     * property evaluation at Eckert reference conditions
	     * and adiabatic wall temperature based on recovery factor.
	     */
	    f = f_flat_plate( Re_L );
	} else {
	    /* 
	     * Default: friction factor determined from 
	     * fully-developed pipe flow correlation.
	     */
	    f = f_darcy_weisbach( Re_D, D ) / lambda;
	}

	if ( L_get_case_id() == DRUMMOND_WITH_M4_NOZZLE && A->sim_time < 0.02 ) {
	    /* After shock-reflection, revert to the pipe-flow friction factor. */
	    f = f_darcy_weisbach( Re_D, D ) / lambda;
	}

#       if DEBUG >= 3
        printf("L_source_vector(): f=%e\n", f);
#       endif

	/* 
	 * Local shear stress, in Pa, computed from the friction factor. 
	 */
	tau0 = -(0.125 * f) * cell->gas->rho * cell->u * abs_u;
	cell->shear_stress = tau0;

	/*
	 * Rate of energy transfer into the cell.
	 * First, shaft work (shear stress).
	 * I think that this term should be zero not the following...
	 * w_dot = tau0 * cell->u * length * mypi * D;
	 */
	w_dot = 0.0;
  
	/*
	 * Convective heat transfer coefficient.
	 */
	St = (f * 0.125) * pow(Prandtl, -0.667);
	h = cell->ref->rho * gmodel->Cp(*(cell->gas)) * abs_u * St;

	/* DJM special */
	if ( A->viscous_effects == VISCOUS_LOSS_PIPE_F_HALF_H ) h *= 0.5;

	/*
	 * Convective heat transfer from the wall into the gas cell.
	 */
	if ( A->adiabatic == 1 ) {
	    q = 0.0;
	} else {
	    q = h * mypi * D * length * (T_wall - T_aw);
	}

	/* 
	 * Compute the energy radiated by the gas cell to the tube wall.
	 */
#       define  RADIATION  0
	if (RADIATION == 1) {
	    Tgas = cell->gas->T[0];
	    alpha = 10.0;   /* 1/metres */
	    sigma = 5.6696e-8;  /* W/(m**2.K**4) */
	    q_rad = 2.0 * length * D * D * mypi *
		alpha * exp(-alpha * 0.5 * D) *
		sigma * Tgas * Tgas * Tgas * Tgas;
	} else {
	    q_rad = 0.0;
	}

	/*
	 * Record the heat flux for the output file.
	 */
	cell->heat_flux = (q - q_rad) / (mypi * D * length);
#       if DEBUG >= 3
	printf( "L_source_vector(): rho=%e, mypi=%e, length=%e, D=%e\n",
		cell->gas->rho, mypi, length, D);
        printf( "L_source_vector(): tau0=%e, q=%e, q_rad=%e\n",
		tau0, q, q_rad);
#       endif

	/*
         * Now, we actually apply the viscous effects to the 
	 * conserved quantities in the cell.
         *
	 * Decide how to remove the mass, momentum and energy
	 * from the cell -- with, or without associated mass loss.
	 *
	 * DEFAULT: Momentum and Energy loss without associated mass loss.
	 * Rate of change of Momentum by shear stress and 
	 * energy loss by convective (and possibly radiative) heat transfer.
	 */
	m_loss   = 0.0;
	mom_loss = -(tau0 * length * mypi * D);
	E_loss   = -(w_dot + q - q_rad);

	if ( A->viscous_effects == CJD_MASS_LOSS_LAMINAR ||
	     A->viscous_effects == CJD_MASS_LOSS_TURBULENT ) {
	    /*
	     * Con Doolan's mass-loss model for shocked gas.
	     * m_loss = mass loss per unit time
	     * This is based on the incompressible flat plate boundary layer
	     * equations and transformed to compresssible via the Howarth 
	     * transformation.
	     */
	    r2r1 = cell->gas->rho / A->init_str->gas->rho;
	    if (cell->L_bar > 0.0 && r2r1 > 1.001) {
		/* Only formulated for strong shocks going forward. */
		if ( A->viscous_effects == CJD_MASS_LOSS_TURBULENT ) {
		    kt = 0.375;
		    m_loss = 0.70 * mypi * D * length * cell->ref->rho * abs_u *
			kt / pow(Re_L, 1.0 / 5.0);
		    mom_loss = 0.889 * m_loss * cell->u;
		    E_loss = gmodel->Cp(*(cell->gas)) * cell->ref->T[0] * m_loss;
		} else { /* Laminar */
		    kl = 5.0;
		    m_loss = 0.333 * mypi * D * length * cell->ref->rho * 
			abs_u * kl / sqrt(Re_L);
		    mom_loss = m_loss * cell->u;
		    E_loss = m_loss * (cell->gas->e[0] + abs_u * abs_u * 0.5);
		}
		if (m_loss < MINIMUM_MASS) { /* Can't remember why. */
		    E_loss   = 0.0;
		    mom_loss = 0.0;
		}
		/* 
		 * *** FIX ME ***
		 * DRB had the following modificaton on the CJD mass loss: 
		 * cell->Q_m -= m_loss * (-tau0 * length * mypi * D / mom_loss);
		 * Is it important?
		 */
	    }
	} else if ( A->viscous_effects == DRB_MASS_LOSS_PIPE_F ||
		    A->viscous_effects == DRB_MASS_LOSS_FLAT_PLATE_F ) {
            /* 
	     * Use default momentum and energy loss (via friction factor).
	     * Obtain mass loss from the momentum loss, noting that
	     * cell velocity is the momentum per unit mass of the core flow. 
	     */ 
	    m_loss = mom_loss / cell->u;
	}

#       if DEBUG >= 3
        printf( "L_source_vector(): viscous effects: m_loss=%e, mom_loss=%e E_loss=%e\n", 
	       m_loss, mom_loss, E_loss);
#       endif
	cell->Q_m   -= m_loss;
	cell->Q_mom -= mom_loss;
	cell->Q_E   -= E_loss;

	/*
	 * Pipe fitting loss is a momentum loss on top of other viscous effects. 
	 */
	F_LOSS = -(cell->K_over_L) * 0.5 * cell->gas->rho * cell->u * abs_u;
	cell->Q_mom += area * length * F_LOSS;

#       if (DEBUG >= 3)
	printf("L_source_vector(): Q_m=%e, Q_mom=%e, Q_E=%e\n", 
	       cell->Q_m, cell->Q_mom, cell->Q_E);
#       endif

    }   /* for (ix = ...  */


#   if DEBUG >= 3
    printf("\nL_source_vector(): End.\n");
#   endif
    return 0;
}   /* end function L_source_vector */


/*-----------------------------------------------------------------*/

int L_axial_heat_flux(struct slug_data *A, double k)
{
    /*
     * Purpose...
     * -------
     * Compute an "axial heat flux" to kill off glitches
     * which seem to appear at the contact surfaces
     * between gas slugs.
     * Note that the form of the flux is quadratic in dT/dx
     * rather than linear.
     * The thermal conductivity, k, is not intended to be accurate. 
     * For air at nominal conditions k = 0.024 SI-units.
     * A value of k = 0.01 seems to work OK for the Sod shock tube.
     */
    int ix;
    double dx, dT;

    for (ix = A->ixmin; ix <= A->ixmax - 1; ++ix) {
        dT = A->Cell[ix + 1].gas->T[0] - A->Cell[ix].gas->T[0];
        dx = A->Cell[ix + 1].xmid - A->Cell[ix].xmid;
        /* Try quadratic, instead of linear... */
        A->Cell[ix].qstar = -k * dT / dx * fabs(dT / dx);
    }   /* end for */

    /* No axial heat transfer to other gas slugs. */
    A->Cell[A->ixmin - 1].qstar = 0.0;
    A->Cell[A->ixmax].qstar = 0.0;

    return 0;
}   /* end function L_axial_heat_flux() */


double min_increment( double y0, double y1, double y2 )
{
    /* Computes the minimum increment to bring y0 into line
       with its neighbour values y1 and y2. */

    double t1, t2, d1, d2, d_min;

    /* Target values */
    t1 = y1;               /* no extrapolation */
    t2 = y1 - (y2 - y1);   /* linear extrapolation */

    /* Increments */
    d1 = t1 - y0;
    d2 = t2 - y0;

    /* Select the one with minimum magnitude. */
    if ( fabs(d1) < fabs(d2) ) {
	d_min = d1;
    } else {
	d_min = d2;
    }

    return d_min;
}


int L_adjust_end_cells(struct slug_data *A)
{
    /*
     * Purpose...
     * -------
     * Adjust the gas properties of the end-of-slug cells
     * to kill off glitches which seem to appear at the 
     * contact surfaces between gas slugs.
     */
    struct L_cell *c, *cn1, *cn2;
    double relax_factor;
    Gas_model *gmodel = get_gas_model_ptr();

    /* 
     * Select a value of the relaxation such that we don't do the 
     * damage all in one go.  
     * Combined with occasional use rather than application every step,
     * this should allow real waves to pass through.
     * However, it doesn't work really well for the start of
     * strong expansions where the pressure and temperature in
     * the first cell *should* drop away really fast.
     */
    relax_factor = 0.1;

    /* 
     * First, do left-end
     */
    c   = &(A->Cell[A->ixmin]);      /* end cell    */
    cn1 = &(A->Cell[A->ixmin + 1]);  /* neighbour 1 */
    cn2 = &(A->Cell[A->ixmin + 2]);  /* neighbour 2 */

    c->gas->T[0] += relax_factor * min_increment(c->gas->T[0], cn1->gas->T[0], cn2->gas->T[0]);
    c->gas->p += relax_factor * min_increment(c->gas->p, cn1->gas->p, cn2->gas->p);
    gmodel->eval_thermo_state_pT(*(c->gas));
    gmodel->eval_transport_coefficients(*(c->gas));
    /* c->u += relax_factor * (cn->u - c->u); */
    L_encode_conserved_for_one_cell( c );

    /* 
     * Second, do right-end
     */
    c   = &(A->Cell[A->ixmax]);      /* end cell    */
    cn1 = &(A->Cell[A->ixmax - 1]);  /* neighbour 1 */
    cn2 = &(A->Cell[A->ixmax - 2]);  /* neighbour 2 */

    c->gas->T[0] += relax_factor * min_increment(c->gas->T[0], cn1->gas->T[0], cn2->gas->T[0]);
    c->gas->p += relax_factor * min_increment(c->gas->p, cn1->gas->p, cn2->gas->p);
    gmodel->eval_thermo_state_pT(*(c->gas));
    gmodel->eval_transport_coefficients(*(c->gas));
    /* c->u += relax_factor * (cn->u - c->u); */
    L_encode_conserved_for_one_cell( c );

    return 0;
}   /* end function L_adjust_end_cells() */

/*-----------------------------------------------------------------*/

int L_apply_rivp(struct slug_data *A)
{
    // Apply the Riemann solver to obtain the pressure and
    // velocity at each interface.
    //
    // Input...
    // A      : pointer to the data structure
    // Output...
    // returns 0 is all is OK, 1 otherwise.
    int ix;
    static struct L_flow_state QL[NDIM], QR[NDIM];
    static double del[NDIM], dplus[NDIM], dminus[NDIM];
    static double rhoL, rhoR;
    static double pstar[NDIM], ustar[NDIM];
    static double onedx[NDIM];
    Gas_model *gmodel = get_gas_model_ptr();
    static double RHOMIN = 1.0e-6;

    // On first encounter, the gas components inside the QL, QR
    // flow state structures will need to be created.
    if ( QL[0].gas == 0 ) {
	// If one is not present, assume that none are present.
	for ( ix = 0; ix < NDIM; ++ix ) {
	    QL[ix].gas = new Gas_data(gmodel);
	    QR[ix].gas = new Gas_data(gmodel);
	}
	printf( "l_apply_rivp(): Setting up workspace for the first time.\n" );
    }

    // Interpolate the cell average values to obtain
    // LEFT and RIGHT states at each interface.
    //
    // The approach taken is to ignore grid distortions and perform
    // one dimensional projection/interpolation in the ix parameter
    // direction.
    // Left cell has index [ix], Right cell has index [ix+1]

    // Always start with first-order interpolation.
    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
	QL[ix].gas->copy_values_from(*(A->Cell[ix].gas));
	QL[ix].u = A->Cell[ix].u;
	QR[ix].gas->copy_values_from(*(A->Cell[ix + 1].gas));
	QR[ix].u = A->Cell[ix + 1].u;
    }

    if ( A->Xorder == 2 ) {
        // Higher-order interoplation of some quantities.
        // Assume 2 ghost points are available.
        for (ix = A->ixmin - 1; ix <= A->ixmax + 2; ++ix) {
            onedx[ix] = 1.0 / (A->Cell[ix].xmid - A->Cell[ix - 1].xmid);
        }
        // Density.
        for (ix = A->ixmin - 1; ix <= A->ixmax + 1; ++ix) {
            dminus[ix] = (A->Cell[ix].gas->rho - A->Cell[ix - 1].gas->rho) * onedx[ix];
            dplus[ix] =
                (A->Cell[ix + 1].gas->rho - A->Cell[ix].gas->rho) * onedx[ix + 1];
            del[ix] = MIN_MOD(dminus[ix], dplus[ix]);
        }
        for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
            rhoL = A->Cell[ix].gas->rho +
                del[ix] * (A->Cell[ix].x - A->Cell[ix].xmid);
            QL[ix].gas->rho = MAXIMUM(rhoL, RHOMIN);
            rhoR = A->Cell[ix + 1].gas->rho -
                del[ix + 1] * (A->Cell[ix + 1].xmid - A->Cell[ix].x);
            QR[ix].gas->rho = MAXIMUM(rhoR, RHOMIN);
        }
        // Axial Velocity
        for (ix = A->ixmin - 1; ix <= A->ixmax + 1; ++ix) {
            dminus[ix] = (A->Cell[ix].u - A->Cell[ix - 1].u) * onedx[ix];
            dplus[ix] = (A->Cell[ix + 1].u - A->Cell[ix].u) * onedx[ix + 1];
            del[ix] = MIN_MOD(dminus[ix], dplus[ix]);
        }
        for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
            QL[ix].u = A->Cell[ix].u +
                del[ix] * (A->Cell[ix].x - A->Cell[ix].xmid);
            QR[ix].u = A->Cell[ix + 1].u -
                del[ix + 1] * (A->Cell[ix + 1].xmid - A->Cell[ix].x);
        }
        // Specific Internal Energy.
	// FIX-ME -- do we need to deal with the other energy modes?
        for (ix = A->ixmin - 1; ix <= A->ixmax + 1; ++ix) {
            dminus[ix] = (A->Cell[ix].gas->e[0] - A->Cell[ix - 1].gas->e[0]) * onedx[ix];
            dplus[ix] = (A->Cell[ix + 1].gas->e[0] - A->Cell[ix].gas->e[0]) * onedx[ix + 1];
            del[ix] = MIN_MOD(dminus[ix], dplus[ix]);
        }
        for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
            QL[ix].gas->e[0] = A->Cell[ix].gas->e[0] +
                del[ix] * (A->Cell[ix].x - A->Cell[ix].xmid);
            QR[ix].gas->e[0] = A->Cell[ix + 1].gas->e[0] -
                del[ix + 1] * (A->Cell[ix + 1].xmid - A->Cell[ix].x);
        }
        // Pressure, Local Speed of Sound and Temperature.
        for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
	    gmodel->eval_thermo_state_rhoe(*(QL[ix].gas));
	    gmodel->eval_thermo_state_rhoe(*(QR[ix].gas));
        }
    } // End of Higher-order interpolation.

    // *****************************
    // * Apply the Riemann solver. *
    // *****************************
    // Apply specified interface velocities if required.
    if (A->set_left_end_ustar == 1) {
        ustar[A->ixmin - 1] = A->left_ustar;
    }
    if (A->set_right_end_ustar == 1) {
        ustar[A->ixmax] = A->right_ustar;
    }
    L_rivp(QL, QR, ustar, pstar, A->ixmin - 1, A->ixmax,
           A->set_left_end_ustar, A->set_right_end_ustar);
    // Save the interface pressures and velocities.
    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
        A->Cell[ix].pface = pstar[ix];
        A->Cell[ix].uface = ustar[ix];
    }
    // Save the interface pressures at the ends in particular.
    A->left_pstar = pstar[A->ixmin - 1];
    A->right_pstar = pstar[A->ixmax];

    return 0;
} // end function L_apply_rivp


/*-----------------------------------------------------------------*/

int L_time_derivatives(struct slug_data *A, int time_level)
{
    /*
     * Purpose...
     * -------
     * Compute the time derivatives for the conserved quantities
     * for each active cell.  These are the spatial (RHS) terms 
     * in the semi-discrete governing equations.
     *
     * 10-Nov-94 : Added a length scale
     * 05-Jun-00 : Axial heat fluxes into energy eqn.
     *
     * Input...
     * -----
     * A          : pointer to THE data structure
     * time_level : specifies where the computed derivatives
     *              are to be stored.
     *
     */
    int ix;
    struct L_cell *C, *Cm1;
    // Interface motions.
    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        C->DxDt[time_level] = C->uface;
    }
    // Cell average properties.
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        Cm1 = &(A->Cell[ix - 1]);
        // Mass.
        C->DmDt[time_level] = C->Q_m;
        // Momentum.
        C->DmomDt[time_level] =
            Cm1->pface * Cm1->area - C->pface * C->area +
            C->gas->p * (C->area - Cm1->area) + C->Q_mom;
	// Energy.
        C->DEDt[time_level] =
            Cm1->pface * Cm1->area * Cm1->uface
            - C->pface * C->area * C->uface
            + Cm1->qstar * Cm1->area - C->qstar * C->area + C->Q_E;
	// Length scale for the application of Mirel's theory.
        C->DLDt[time_level] = fabs(C->u);
    } // end for (ix = ...
    return SUCCESS;
}   /* end function L_time_derivatives  */

//-----------------------------------------------------------------

int L_record_slug_state(struct slug_data *A)
{
    /*
     * Purpose...
     * -------
     * Record the slug state before attempting a time step.
     *
     * 26-April-98 : first introduced
     *
     * Input...
     * -----
     * A      : pointer to the slug_data structure
     *
     */
    int ix;
    struct L_cell *C;

#   if (DEBUG >= 2)
    printf("Record slug state...\n");
#   endif

    /*
     * Interfaces.
     */
    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        C->x_old = C->x;
    }   /* end for */

    /*
     * Conserved quantities within the cell.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        C->mass_old = C->mass;
        C->moment_old = C->moment;
        C->Energy_old = C->Energy;
        C->L_bar_old = C->L_bar;
    }   /* end for */

    return 0;
}   /* end function L_record_slug_state */

/*-----------------------------------------------------------------*/

int L_restore_slug_state(struct slug_data *A)
{
    /*
     * Purpose...
     * -------
     * Restore the cell state to that before 
     * the attempted time step.
     *
     * 26-April-98 : first introduced
     *
     * Input...
     * -----
     * A      : pointer to the slug_data structure
     *
     */
    int ix;
    struct L_cell *C;

#   if (DEBUG >= 2)
    printf("Restore slug state...\n");
#   endif

    /*
     * Interfaces.
     */
    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        C->x = C->x_old;
    }   /* end for */

    /*
     * Conserved quantities within the cell.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        C->mass = C->mass_old;
        C->moment = C->moment_old;
        C->Energy = C->Energy_old;
        C->L_bar = C->L_bar_old;
    }   /* end for */

    return 0;
}   /* end function L_restore_slug_state */

/*-----------------------------------------------------------------*/

int L_predictor_step(struct slug_data *A)
{
    /*
     * Purpose...
     * -------
     * Use the time derivatives to advance the conserved quantities
     * forward by time step dt.
     *
     * 10-Nov-94 : Added a length scale
     *             and lower bound on cell mass
     *
     * Input...
     * -----
     * A      : pointer to the slug_data structure
     *
     */
    int ix;
    double mass, dm, dt;
    struct L_cell *C;

#   if (DEBUG >= 2)
    printf("Predictor step update...\n");
#   endif

    dt = A->dt;

    /*
     * Interfaces.
     */
    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        C->x = C->x_old + dt * C->DxDt[0];
    }

    /*
     * Conserved quantities within the cell.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);

        mass = C->mass_old;
        dm = A->dt * C->DmDt[0];
        if (mass + dm > MINIMUM_MASS)
            C->mass = mass + dm;
	else
	    C->mass = mass;
        C->moment = C->moment_old + dt * C->DmomDt[0];
        C->Energy = C->Energy_old + dt * C->DEDt[0];

        C->L_bar = C->L_bar_old + dt * C->DLDt[0];
    }   /* end for */

    return 0;
}   /* end function L_predictor_step */

/*-----------------------------------------------------------------*/

int L_corrector_step(struct slug_data *A)
{
    /*
     * Purpose...
     * -------
     * Use the time derivatives to advance the conserved quantities
     * forward by time step dt.
     *
     * 10-Nov-94 : Added length scale.
     *             lower bound on cell mass
     * 26-Apr-98 : clean up pointers
     *
     * Input...
     * -----
     * A      : pointer to the slug_data structure
     *
     */
    int ix;
    double mass, dm, dt;
    struct L_cell *C;

#   if (DEBUG >= 2)
    printf("Corrector step update...\n");
#   endif

    dt = A->dt;

    /*
     * Interfaces.
     */
    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        C->x = C->x_old + dt * 0.5 * (C->DxDt[1] + C->DxDt[0]);
    }   /* end for */

    /*
     * Conserved quantities within the cell.
     */
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);

        mass = A->Cell[ix].mass_old;
        dm = dt * 0.5 * (C->DmDt[1] + C->DmDt[0]);
        if (mass + dm > MINIMUM_MASS)
            C->mass = mass + dm;
	else
	    C->mass = mass;
        C->moment = C->moment_old + dt * 0.5 * (C->DmomDt[1] + C->DmomDt[0]);
        C->Energy = C->Energy_old + dt * 0.5 * (C->DEDt[1] + C->DEDt[0]);
        /* length scale */
        C->L_bar = C->L_bar_old + dt * 0.5 * (C->DLDt[1] + C->DLDt[0]);
    }   /* end for */

    return 0;
}   /* end function L_corrector_step */

/*-----------------------------------------------------------------*/

int L_check_cells(struct slug_data *A, int js)
{
    /*
     * Purpose...
     * -------
     * Check cells within a given slug for invalid data.
     *
     * Input...
     * -----
     * A      : pointer to the slug_data structure
     * js     : index for the particular slug
     *
     * Output...
     * ------
     * L_check_cells() returns the number of bad cells.
     *
     */
    int ix, number_bad_cells, this_cell_bad;
    struct L_cell *C, *Cm1;
    double dx;

#   if (DEBUG >= 2)
    printf("Check for bad cells...\n");
#   endif

    number_bad_cells = 0;
    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        C = &(A->Cell[ix]);
        Cm1 = &(A->Cell[ix - 1]);

        this_cell_bad = 0;

#       if DEBUG >= 3
        printf("Slug %d Cell %d: uL=%e, u=%e, uR=%e\n",
               js, ix - A->ixmin, Cm1->uface, C->u, C->uface);
        printf("               : pL=%e, p=%e, pR=%e\n",
               Cm1->pface, C->gas->p, C->pface);
#       endif

        dx = C->x - Cm1->x;
        if (dx <= 0.0)
            this_cell_bad = 1;

        if (C->mass <= 0.0) this_cell_bad = 1;
        if (C->gas->T[0] <= 1.0) this_cell_bad = 1;
        if (C->gas->p <= 1.0e-6) this_cell_bad = 1;

        if (this_cell_bad == 1) {
            ++number_bad_cells;
            printf("Bad Cell: js=%d, cell=%d\n", js, ix - A->ixmin);
            printf("        : dx=%e, xL=%e, xR=%e\n", dx, Cm1->x, C->x);
            printf("        : mass=%e, moment=%e, E=%e\n",
                   C->mass, C->moment, C->Energy);
            printf("        : mass_old=%e, moment_old=%e, E_old=%e\n",
                   C->mass_old, C->moment_old, C->Energy_old);
            printf("        : uLface=%e, u=%e, uRface=%e\n",
                   Cm1->uface, C->u, C->uface);
            printf("        : pLface=%e, p=%e, pRface=%e\n",
                   Cm1->pface, C->gas->p, C->pface);
            printf("        : ALface=%e, volume=%e, ARface=%e\n",
                   Cm1->area, C->volume, C->area);
        }   /* end if */
    }   /* end for */

    return number_bad_cells;
}   /* end function L_check_bad_cells() */

/*--------------------------------------------------*/

int L_check_cfl(struct slug_data *A)
{
    /* 
     * Purpose...
     * -------
     * Compute the local time step limit for each cell.
     * The overall time step is limited by the worst-case cell.
     *
     * Input...
     * -----
     * A        : pointer to the slug_data structure
     *
     * Some Definitions...
     * ----------------
     * dt       : global time step for the block
     * cfl_target : desired CFL number
     * cfl_min  : approximate minimum CFL number in the block
     * cfl_max  : approximate maximum CFL number in the block
     * dt_allow : allowable time step (i.e. the maximum dt that
     *            satisfies the CFL target)
     *
     * 07-dec-00 : check the magnitudes of the interface velocities
     *             rather than just the cell bulk velocity. 
     */
    int ix;
    double signal_time, dt_local, cfl_local;
    double uL_mag, uR_mag, u_mag, a_local, signal_speed, length;

#   if (DEBUG >= 2)
    printf("Checking CFL number...\n");
#   endif

    A->cfl_min = 1.0e6; /* something outrageously large */
    A->cfl_max = 0.0;
    A->dt_allow = 1.0e6;

    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
        /* With the Lagrangian formulation,
         * signal speed is essentially the speed of sound. */
        a_local = A->Cell[ix].gas->a;
        signal_speed = a_local;

#       if 0
        /* Select the largest signal speed from 
         * the velocities, also. */
        u_mag = fabs(A->Cell[ix].u);
        uL_mag = fabs(A->Cell[ix - 1].uface);
        uR_mag = fabs(A->Cell[ix].uface);
        if (signal_speed < uL_mag)
            signal_speed = uL_mag;
        if (signal_speed < uR_mag)
            signal_speed = uR_mag;
        if (signal_speed < u_mag)
            signal_speed = u_mag;
#       else
	UNUSED_VARIABLE(u_mag);
	UNUSED_VARIABLE(uL_mag);
	UNUSED_VARIABLE(uR_mag);
#       endif

        /* 
         * Check the INVISCID time step limit.
         */
        length = A->Cell[ix].x - A->Cell[ix - 1].x;
        signal_time = length / signal_speed;

        /*
         * Current (Local) CFL number
         */
        cfl_local = fabs(A->dt / signal_time);
        /* 
         * Recommend a time.
         */
        dt_local = A->cfl_target * signal_time;

        /*
         * Search for the worst case.
         */
        if (cfl_local < A->cfl_min) A->cfl_min = cfl_local;
        if (cfl_local > A->cfl_max) A->cfl_max = cfl_local;
        if (dt_local < A->dt_allow) A->dt_allow = dt_local;

        /*
         * Some debug for problem situations ...
         */
        if (cfl_local > 10.0) {
            printf("\n-----------------\n");
            printf("CFL = %e , cell[%d]\n", cfl_local, ix);
            printf("rho = %e, u = %e, a = %e, length = %e\n",
                   A->Cell[ix].gas->rho, A->Cell[ix].u, 
		   A->Cell[ix].gas->a, length);
        }   /* end if */
    }   /* end for ix loop */

    if (A->cfl_max > 0.9) {
        printf("WARNING: large CFL number, cfl_max = %e\n", A->cfl_max);
    }

    return 0;
}   /* end function L_check_cfl */

/*=================================================================*/

int P_time_derivatives(struct piston_data *B, int time_level, double sim_time)
{
    /*
     * Purpose...
     * -------
     * Compute the time derivatives for the piston dynamics
     *
     * Input...
     * -----
     * B          : pointer to the piston_data structure
     * time_level : specifies where the computed derivatives
     *              are to be stored.
     * sim_time   : current simulation time
     *
     */
    double pressure_force, friction_force;
    /* The tolerance below which the velocity is assumed zero. */
#   define V_TOL 1.0e-6

    /*
     * If we have a restrained piston, check for release.
     */
    if (B->is_restrain == 1) {
        if (B->Pb > B->p_restrain)
            B->is_restrain = 0;
    }   /* end if */

    /*
     * If the piston is still restrained, we are going nowhere.
     */
    if (B->is_restrain == 1) {
        /* Set speed and acceleration to zero. */
        B->V = 0.0;
        B->DxDt[time_level] = 0.0;
        B->DVDt[time_level] = 0.0;
        return 0;
    }   /* end if */


    /*
     * Check for buffer strike only while we're moving forward.
     */
    if (B->x > B->x_buffer && B->V > V_TOL && B->on_buffer == 0) {
        B->hit_buffer_count += 1;
        B->on_buffer = 1;

        /* Record striking speed. */
        B->V_buffer = B->V;
        printf("Buffer strike: speed = %e\n", B->V_buffer);
    }   /* end if */

    /*
     * If we were on the buffer and the velocity is now negative,
     * we must have moved off the buffer.
     */
    if (B->on_buffer == 1 && B->V < 0.0) {
        B->on_buffer = 0;
    }   /* end if */

    /*
     * Consequence of being on the buffer is that we have zero speed.
     * There may still be pressure and friction forces to consider.
     */
    if (B->on_buffer == 1) {
        /* Set speed to zero for this step. */
        B->V = 0.0;
    }


    /*
     * Check for the application of brakes.
     */
    if ((B->with_brakes == 1) && (B->V < 0.0))
        B->brakes_on = 1;

    /*
     * If the brakes are on, we are going nowhere.
     */
    if (B->brakes_on == 1) {
        /* Set speed and acceleration to zero. */
        B->V = 0.0;
        B->DxDt[time_level] = 0.0;
        B->DVDt[time_level] = 0.0;
        return 0;
    }   /* end if */

    /*
     * ***********************************************
     * Now continue on with a free-piston calculation.
     * ***********************************************
     *
     * The (signed) pressure force.
     */
    pressure_force = B->area * (B->Pb - B->Pf);
    /*
     * The magnitude of the friction force.
     */
    friction_force = B->front_seal_f * B->front_seal_area * B->Pf +
	B->back_seal_f * B->back_seal_area * B->Pb;


    /* 
     * Update the state vector for the piston dynamics.
     *
     * This is the velocity.
     */
    B->DxDt[time_level] = B->V;


    /* 
     * Now, do the acceleration.
     */
    if (B->V > V_TOL) {
	/* 
	 * Moving forward, apply full friction in reverse 
	 */
	B->DVDt[time_level] = (pressure_force - friction_force) / B->mass;

    } else if (B->V < -(V_TOL)) {
	/* 
	 * Moving backward, apply full friction forward 
	 */
	B->DVDt[time_level] = (pressure_force + friction_force) / B->mass;

    } else {
	/*
	 * We are effectively stationary.
	 *
	 * If the pressure force is larger than the friction
	 * force then apply the difference else remain stationary.
	 */
	if (fabs(pressure_force) > friction_force) {
	    /* 
	     * Pressure force dominates, apply the remainder 
	     */
	    if (pressure_force > 0.0) {
		/* Accelerate forward. */
		B->DVDt[time_level] =
		    (pressure_force - friction_force) / B->mass;
	    } else {
		/* Accelerate backwards */
		B->DVDt[time_level] =
		    (pressure_force + friction_force) / B->mass;
	    } /* end if */

	} else {
	    /* 
	     * Friction force dominates; let's remain stationary. 
	     */
	    B->DVDt[time_level] = 0.0;
	} /* end if */
    } /* end if ... stationary */

    return SUCCESS;
} /* end function P_time_derivatives */

/*-----------------------------------------------------------------*/

int P_record_piston_state(struct piston_data *B)
{
    /*
     * Purpose...
     * -------
     * Record the piston position and velocity before
     * attempting a time-step.
     *
     * Input...
     * -----
     * B      : pointer to the piston_data structure
     *
     */
#   if (DEBUG >= 2)
    printf("Piston record state...\n");
#   endif

    B->x_old = B->x;
    B->V_old = B->V;
    B->brakes_on_old = B->brakes_on;
    B->hit_buffer_count_old = B->hit_buffer_count;
    B->on_buffer_old = B->on_buffer;
    B->is_restrain_old = B->is_restrain;
    B->mass_old = B->mass;

    return 0;
}   /* end function P_record_piston_state */

/*-----------------------------------------------------------------*/

int P_restore_piston_state(struct piston_data *B)
{
    /*
     * Purpose...
     * -------
     * Restore the piston position and velocity to that 
     * which existed before attempting the time-step.
     *
     * Input...
     * -----
     * B      : pointer to the piston_data structure
     *
     */
#   if (DEBUG >= 2)
    printf("Piston restore state...\n");
#   endif

    B->x = B->x_old;
    B->V = B->V_old;
    B->brakes_on = B->brakes_on_old;
    B->hit_buffer_count = B->hit_buffer_count_old;
    B->on_buffer = B->on_buffer_old;
    B->is_restrain = B->is_restrain_old;
    B->mass = B->mass_old;

    return 0;
}   /* end function P_restore_piston_state */

/*-----------------------------------------------------------------*/

int P_predictor_step(struct piston_data *B)
{
    /*
     * Purpose...
     * -------
     * Use the time derivatives to advance the piston dynamics
     * forward by time step dt.
     *
     * Input...
     * -----
     * B      : pointer to the piston_data structure
     *
     */
#   if (DEBUG >= 2)
    printf("Piston Predictor step update...\n");
#   endif

    if (B->is_restrain || B->brakes_on) {
        B->x = B->x_old;
        B->V = 0.0;
    } else if (B->on_buffer == 1) {
        B->x = B->x_old;
        B->V = 0.0 + B->dt * B->DVDt[0];
    } else {
        B->x = B->x_old + B->dt * B->DxDt[0];
        B->V = B->V_old + B->dt * B->DVDt[0];
    }   /* end if */
    
    /*
     * Apply mass decay for flagged and unrestrained pistons only.
     */
    
    if (B->is_restrain == 0 && B->apply_decay == 1) {
        /* Limit the mass decay to a user-defined level */
        if (B->mass > B->mass_limit) {
            B->mass = B->mass_old * ( 1.0 - B->f_decay * B->dt );
        }
        else {
            B->mass = B->mass_limit;
        }
    }

    return 0;
}   /* end function P_predictor_step */

/*-----------------------------------------------------------------*/

int P_corrector_step(struct piston_data *B)
{
    /*
     * Purpose...
     * -------
     * Use the time derivatives to advance the piston dynamics
     * forward by time step dt.
     *
     * Input...
     * -----
     * B      : pointer to the piston_data structure
     *
     */
#   if (DEBUG >= 2)
    printf("Piston Corrector step update...\n");
#   endif

    if (B->is_restrain || B->brakes_on) {
        B->x = B->x_old;
        B->V = 0.0;
    } else if (B->on_buffer == 1) {
        B->x = B->x_old;
        B->V = 0.0 + B->dt * 0.5 * (B->DVDt[1] + B->DVDt[0]);
    } else {
        B->x = B->x_old + B->dt * 0.5 * (B->DxDt[1] + B->DxDt[0]);
        B->V = B->V_old + B->dt * 0.5 * (B->DVDt[1] + B->DVDt[0]);
    }   /* end if */
    
    return 0;
}   /* end function P_corrector_step */


/*================= end of l_tstep.c ==================*/
