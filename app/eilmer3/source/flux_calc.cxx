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
#include "../../../lib/util/source/useful.h"
#include "../../../lib/geometry2/source/geom.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"


// Local flow_state object for temporarily holding interface state
// while computing flux at each interface.
FlowState *IFace_flow_state = NULL;


std::map<std::string,flux_calc_t> available_calculators = {
    {"0",FLUX_RIEMANN}, {"riemann",FLUX_RIEMANN}, {"Riemann",FLUX_RIEMANN},
    {"1",FLUX_AUSM}, {"ausm",FLUX_AUSM}, {"AUSM",FLUX_AUSM},
    {"2",FLUX_EFM}, {"efm",FLUX_EFM}, {"EFM",FLUX_EFM},
    {"3",FLUX_AUSMDV}, {"ausmdv",FLUX_AUSMDV}, {"AUSMDV",FLUX_AUSMDV},
    {"4",FLUX_ADAPTIVE}, {"adaptive",FLUX_ADAPTIVE}, {"Adaptive",FLUX_ADAPTIVE},
    {"5",FLUX_AUSM_PLUS_UP}, {"ausm_plus_up",FLUX_AUSM_PLUS_UP}, {"AUSM_plus_up",FLUX_AUSM_PLUS_UP},
    {"6",FLUX_HLLE}, {"hlle",FLUX_HLLE}, {"HLLE",FLUX_HLLE},
};

std::map<flux_calc_t,std::string> calculator_names = {
    {FLUX_RIEMANN,"riemann"},
    {FLUX_AUSM,"ausm"},
    {FLUX_EFM,"efm"},
    {FLUX_AUSMDV,"ausmdv"},
    {FLUX_ADAPTIVE,"adaptive"},
    {FLUX_AUSM_PLUS_UP,"ausm_plus_up"},
    {FLUX_HLLE,"hlle"}
};

flux_calc_t flux_calculator = FLUX_RIEMANN;

flux_calc_t set_flux_calculator(std::string name)
{
    flux_calculator = available_calculators[name];
    return flux_calculator;
}

flux_calc_t get_flux_calculator()
{
    return flux_calculator;
}

std::string get_flux_calculator_name(flux_calc_t calc)
{
    return calculator_names[calc];
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
    double WSL, WSR;
    ConservedQuantities &F = *(IFace.F);
    Gas_model *gmodel = get_gas_model_ptr();
    if( IFace_flow_state == NULL ) {
	IFace_flow_state = new FlowState(gmodel);
    }

    // IFace.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    // Lft.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    // Rght.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);

    // Transform to interface frame of reference.
    Lft.vel -= IFace.vel;
    Rght.vel -= IFace.vel;
    IFace.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    Lft.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    Rght.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);

    // also transform the magnetic field
    if (get_mhd_flag() == 1) {
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
        throw runtime_error("Invalid flux calculator.");
    } // end switch

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

int set_flux_vector_in_local_frame(ConservedQuantities &F, const FlowState &fs)
{
    double rho = fs.gas->rho;
    double un = fs.vel.x;
    double vt1 = fs.vel.y;
    double vt2 = fs.vel.z;
    double p = fs.gas->p;
    double e = fs.gas->e[0];
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
    ConservedQuantities &F = *(IFace.F);
    double vx = fs.vel.x; double vy = fs.vel.y; double vz = fs.vel.z; // to restore fs at end
    // Transform to interface frame of reference.
    fs.vel -= IFace.vel; // Beware: fs.vel is changed here and restored below.
    IFace.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    fs.vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    // also transform the magnetic field
    if (get_mhd_flag() == 1) {
	fs.B.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    }
    set_flux_vector_in_local_frame(*(IFace.F), fs);
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	double x = IFace.pos.x;
	double y = IFace.pos.y;
	double rsq = x*x + y*y;
	// The conserved quantity is rothalpy,
	// so we need to take -(u**2)/2 off the total energy Shock.
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
    if (get_mhd_flag() == 1) {
	F.B.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    }
    fs.vel.x = vx; fs.vel.y = vy; fs.vel.z = vz; // restore fs.vel
    return SUCCESS;
} // end compute_boundary_flux()
