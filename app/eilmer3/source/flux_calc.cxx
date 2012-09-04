/** \file flux_calc.cxx
 * \ingroup eilmer3
 * \brief Generic flux calculation function for Elmer3, etc.
 *
 * \author PA Jacobs
 *
 * \version 05-Aug-04 : Extracted from mb_cns/source/cns_invs.c.
 * ]version Feb-Jul-2008  : Elmer3 port
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/geometry2/source/geom.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"

/*--------------------------------------------------------------*/

int flux_calculator;    /* 0 == FLUX_RIEMANN */
                        /* 1 == FLUX_AUSM    */
                        /* 2 == FLUX_EFM     */
                        /* 3 == FLUX_AUSMDV  */
                        /* 4 == FLUX_ADAPTIVE */
                        /* 5 == FLUX_AUSM_PLUS_UP */
                        /* 6 == FLUX_HLLE */
/* Local flow_state structures for temporarily holding interface state */
FlowState *IFace_flow_state = NULL;


/*--------------------------------------------------------------*/

/** \brief Sets a global variable to select the flux calculator. */
int set_flux_calculator(int iflux)
{
    flux_calculator = iflux;
    if (flux_calculator == FLUX_RIEMANN)
        if ( get_verbose_flag() ) printf("Fluxes calculated via Riemann solver\n");
    else if (flux_calculator == FLUX_AUSM)
        if ( get_verbose_flag() ) printf("Fluxes calculated via AUSM\n");
    else if (flux_calculator == FLUX_EFM)
        if ( get_verbose_flag() ) printf("Fluxes calculated via EFM\n");
    else if (flux_calculator == FLUX_AUSMDV)
        if ( get_verbose_flag() ) printf("Fluxes calculated via AUSMDV\n");
    else if (flux_calculator == FLUX_ADAPTIVE)
        if ( get_verbose_flag() )  printf("Fluxes calculated via ADAPTIVE\n");
    else if (flux_calculator == FLUX_AUSM_PLUS_UP)
        if ( get_verbose_flag() ) printf("Fluxes calculated via AUSM_PLUS_UP\n");
    else if (flux_calculator == FLUX_HLLE)
        if ( get_verbose_flag() ) printf("Fluxes calculated via HLLE\n");
    else {
        printf("Invalid flux calculator: %d. Riemann solver assumed\n",
	       flux_calculator);
        flux_calculator = FLUX_RIEMANN;
    }
    return 0;
}


/** \brief Returns the index of the selected flux calculator. */
int get_flux_calculator(void)
{
    return flux_calculator;
}

/*----------------------------------------------------------------------------*/

/** \brief Compute the inviscid fluxes (in 2D) across the cell interfaces.
 *
 * This is the top-level function that calls the previously selected
 * flux calculator.
 * Much of the detailed work is delegated.
 *
 * \param Lft    : reference to the LEFT flow state
 * \param Rght   : reference to the RIGHT flow state
 * \param IFace  : reference to the interface where the fluxes are to be stored
 */
int compute_interface_flux(FlowState &Lft, FlowState &Rght, FV_Interface &IFace, double omegaz)
{
    double rho, e, p, ke;
    double un, vt1, vt2;
    double WSL, WSR;
    ConservedQuantities &F = *(IFace.F);
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = get_gas_model_ptr()->get_number_of_species();
    int nmodes = get_gas_model_ptr()->get_number_of_modes();
    if( IFace_flow_state == NULL ) {
	IFace_flow_state = new FlowState(gmodel);
    }

    Lft.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    Rght.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);

    // also transform the magnetic field
    if (get_mhd_flag() == 1) {
        Lft.B.transform_to_local(IFace.n, IFace.t1, IFace.t2);
        Rght.B.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    }

    /*
     * Compute the fluxes in the local frame of the interface.
     */
    if (get_flux_calculator() == FLUX_AUSM ||
        get_flux_calculator() == FLUX_RIEMANN) {
        if (get_flux_calculator() == FLUX_RIEMANN)
            rivp(Lft, Rght, *IFace_flow_state, WSL, WSR);
        else
            ausm(Lft, Rght, *IFace_flow_state, WSL, WSR);

	rho = IFace_flow_state->gas->rho;
	un = IFace_flow_state->vel.x;
	vt1 = IFace_flow_state->vel.y;
	vt2 = IFace_flow_state->vel.z;
	p = IFace_flow_state->gas->p;
	e = IFace_flow_state->gas->e[0];
	ke = 0.5 * (un*un + vt1*vt1 + vt2*vt2);
	/* Kinetic energy per unit volume. */

	/* Mass flux (mass / unit time / unit area) */
	F.mass = rho * un;
	/* Flux of normal momentum */
	F.momentum.x = rho * un * un + p;
	/* Flux of tangential momentum */
	F.momentum.y = rho * un * vt1;
	F.momentum.z = rho * un * vt2;
	/* Flux of Total Energy */
	F.total_energy = rho * un * (e + ke) + p * un;
	F.tke = rho * un * IFace_flow_state->tke;  // turbulence kinetic energy
	F.omega = rho * un * IFace_flow_state->omega;  // pseudo vorticity
	/* Species mass flux */
	for ( int isp = 0; isp < nsp; ++isp ) {
	    F.massf[isp] = F.mass * IFace_flow_state->gas->massf[isp];
	}
	/* Individual energies. */
	// NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = F.mass * IFace_flow_state->gas->e[imode];
	}
    } else if (get_flux_calculator() == FLUX_EFM) {
        efmflx(Lft, Rght, IFace);
    } else if (get_flux_calculator() == FLUX_AUSMDV) {
        ausmdv(Lft, Rght, IFace);
    } else if (get_flux_calculator() == FLUX_ADAPTIVE) {
        adaptive_flux(Lft, Rght, IFace);
    } else if (get_flux_calculator() == FLUX_AUSM_PLUS_UP) {
        ausm_plus_up(Lft, Rght, IFace);
    } else if (get_flux_calculator() == FLUX_HLLE) {
        hlle(Lft, Rght, IFace);
    } else {
        printf("Invalid flux calculator.\n");
        exit(VALUE_ERROR);
    }

    if ( omegaz != 0.0 ) {
	// Rotating frame.
	double x = IFace.pos.x;
	double y = IFace.pos.y;
	double rsq = x*x + y*y;
	// The conserved quantity is rothalpy,
	// so we need to take -(u**2)/2 off the total energy flux.
	// Note that rotating frame velocity u = omegaz * r.
	F.total_energy -= F.mass * 0.5*omegaz*omegaz*rsq;
    }

#   if DEBUG_FLUX >= 1
    printf("Inviscid Fluxes in local frame\n");
    F.print();
#   endif

    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	
	// also transform the magnetic field
    if (get_mhd_flag() == 1) {
      F.B.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    }
	
#   if DEBUG_FLUX >= 1
    printf("Interface fluxes\n");
    printf("xyz_mom.x=%e, \nxyz_mom.y=%e, xyz_mom.z=%e\n", F.momentum.x, F.momentum.y, F.momentum.z);
    if (get_mhd_flag() == 1) {
      printf("xyz_B.x=%e, \nxyz_B.y=%e, xyz_B.z=%e\n", F.B.x, F.B.y, F.B.z);
    }
#   endif

    return SUCCESS;
} // end of compute_interface_flux()

