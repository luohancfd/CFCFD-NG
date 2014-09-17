/** \file flux_calc.cxx
 * \ingroup eilmer3
 * \brief Generic flux calculation function for Elmer3, etc.
 *
 * \author PA Jacobs
 *
 * \version 05-Aug-04 : Extracted from mb_cns/source/cns_invs.c.
 * ]version Feb-Jul-2008  : Elmer3 port
 */

#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/geometry2/source/geom.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"


// Local flow_state object for temporarily holding interface state
// while computing flux at each interface.
FlowState *IFace_flow_state = NULL;

flux_calc_t flux_calculator = FLUX_RIEMANN;

flux_calc_t set_flux_calculator(flux_calc_t my_flux_calculator)
{
    flux_calculator = my_flux_calculator;
    return flux_calculator;
}

flux_calc_t get_flux_calculator()
{
    return flux_calculator;
}

std::string get_flux_calculator_name(flux_calc_t calc)
{
    switch ( calc ) {
    case FLUX_RIEMANN: return "riemann";
    case FLUX_AUSM: return "ausm";
    case FLUX_EFM: return "efm";
    case FLUX_AUSMDV: return "ausmdv";
    case FLUX_ADAPTIVE: return "adaptive";
    case FLUX_AUSM_PLUS_UP: return "ausm_plus_up";
    case FLUX_HLLE: return "hlle";
    default: return "none";
    }
} // end get_flux_calculator_name()

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
    global_data &G = *get_global_data_ptr();
    double WSL, WSR;
    ConservedQuantities &F = *(IFace.F);
    Gas_model *gmodel = get_gas_model_ptr();
    if( IFace_flow_state == NULL ) {
	IFace_flow_state = new FlowState(gmodel);
    }

    // Transform to interface frame of reference.
    Lft.vel.x -= IFace.vel.x;  Lft.vel.y -= IFace.vel.y;  Lft.vel.z -= IFace.vel.z;
    Rght.vel.x -= IFace.vel.x; Rght.vel.y -= IFace.vel.y; Rght.vel.z -= IFace.vel.z;
    IFace.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    Lft.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    Rght.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);

    // also transform the magnetic field
    if ( G.MHD ) {
	Lft.B.transform_to_local(IFace.n, IFace.t1, IFace.t2);
        Rght.B.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    }

    // Compute the fluxes in the local frame of the interface.
    switch ( get_flux_calculator() ) {
    case FLUX_AUSM:
    case FLUX_RIEMANN:
        if (get_flux_calculator() == FLUX_RIEMANN) {
            rivp(Lft, Rght, *IFace_flow_state, WSL, WSR);
        } else {
            ausm(Lft, Rght, *IFace_flow_state, WSL, WSR);
	}
	set_flux_vector_in_local_frame(*(IFace.F), *IFace_flow_state);
	break;
    case FLUX_EFM:
        efmflx(Lft, Rght, IFace);
	break;
    case FLUX_AUSMDV:
        ausmdv(Lft, Rght, IFace);
	break;
    case FLUX_ADAPTIVE:
        adaptive_flux(Lft, Rght, IFace);
	break;
    case FLUX_AUSM_PLUS_UP:
        ausm_plus_up(Lft, Rght, IFace);
	break;
    case FLUX_HLLE:
        hlle(Lft, Rght, IFace);
	break;
    default:
        throw std::runtime_error("Invalid flux calculator.");
    } // end switch

    // perform divergence cleaning of the magnetic field
    // implemented by modifying the flux of B according to the method of Dedner et al.
    if ( G.MHD ) {

	F.divB = Lft.B.x +  0.5 * (Rght.B.x - Lft.B.x) - 0.5 * (Rght.psi - Lft.psi) / G.c_h;

	F.B.x += Lft.psi +  0.5 * (Rght.psi - Lft.psi) - 0.5 * G.c_h * (Rght.B.x - Lft.B.x);
	F.psi  = G.c_h * G.c_h * F.divB;
	
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
    
    constexpr bool do_print_fluxes = false; // just for debug
    if ( do_print_fluxes ) {
	printf("Inviscid Fluxes in local frame\n");
	F.print();
    }

    // Transform fluxes back from interface frame of reference to local frame of reference.
    /* Flux of Total Energy */
    F.total_energy += 0.5 * F.mass * pow(vabs(IFace.vel),2) + dot(F.momentum, IFace.vel);
    /* Flux of momentum */
    F.momentum += F.mass * IFace.vel;
 
    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    // also transform the interface velocities
    IFace.vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    // also transform the magnetic field
    if ( G.MHD ) {
	F.B.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    }
	
    if ( do_print_fluxes ) {
	printf("Interface fluxes\n");
	printf("xyz_mom.x=%e, \nxyz_mom.y=%e, xyz_mom.z=%e\n", 
	       F.momentum.x, F.momentum.y, F.momentum.z);
	if ( G.MHD ) {
	    printf("xyz_B.x=%e, \nxyz_B.y=%e, xyz_B.z=%e\n", 
		   F.B.x, F.B.y, F.B.z);
	}
    }
    
    return SUCCESS;
} // end of compute_interface_flux()

int set_flux_vector_in_local_frame(ConservedQuantities &F, const FlowState &fs)
{
    double rho = fs.gas->rho;
    double un = fs.vel.x;
    double vt1 = fs.vel.y;
    double vt2 = fs.vel.z;
    double p = fs.gas->p;
    double e = accumulate(fs.gas->e.begin(), fs.gas->e.end(), 0.0);
    double ke = 0.5 * (un*un + vt1*vt1 + vt2*vt2); // Kinetic energy per unit volume.
    
    // Mass flux (mass / unit time / unit area)
    F.mass = rho * un; // The mass flux is relative to the moving interface.
    // Flux of normal momentum
    F.momentum.x = F.mass * un + p;
    // Flux of tangential momentum
    F.momentum.y = F.mass * vt1;
    F.momentum.z = F.mass * vt2;
    // Flux of Total Energy
    F.total_energy = F.mass * (e + ke) + p * un;
    F.tke = F.mass * fs.tke;  // turbulence kinetic energy
    F.omega = F.mass * fs.omega;  // pseudo vorticity
    // Species mass flux
    for ( size_t isp = 0; isp < F.massf.size(); ++isp ) {
	F.massf[isp] = F.mass * fs.gas->massf[isp];
    }
    // Individual energies.
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    for ( size_t imode = 1; imode < F.energies.size(); ++imode ) {
	F.energies[imode] = F.mass * fs.gas->e[imode];
    }
    return SUCCESS;
}

int set_flux_vector_in_global_frame(FV_Interface &IFace, FlowState &fs, double omegaz)
{
    global_data &G = *get_global_data_ptr();
    ConservedQuantities &F = *(IFace.F);
    double vx = fs.vel.x; double vy = fs.vel.y; double vz = fs.vel.z; // to restore fs at end
    // Transform to interface frame of reference.
    fs.vel -= IFace.vel; // Beware: fs.vel is changed here and restored below.
    IFace.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    fs.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    // also transform the magnetic field
    if ( G.MHD ) {
	fs.B.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    }
    set_flux_vector_in_local_frame(*(IFace.F), fs);
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

    // Transform fluxes back from interface frame of reference to local frame of reference.
    /* Flux of Total Energy */
    F.total_energy += 0.5 * F.mass * pow(vabs(IFace.vel),2) + dot(F.momentum, IFace.vel);
    /* Flux of momentum */
    F.momentum += F.mass * IFace.vel;

    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    // also transform the interface velocities
    IFace.vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);	  
    // also transform the magnetic field
    if ( G.MHD ) {
	F.B.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    }
    fs.vel.x = vx; fs.vel.y = vy; fs.vel.z = vz; // restore fs.vel
    return SUCCESS;
} // end compute_boundary_flux()
