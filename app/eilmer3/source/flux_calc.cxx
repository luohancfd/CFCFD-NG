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
    case FLUX_HLLC: return "hllc";
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
    // IFace.ivel is the moving velocity of interface
    Lft.vel.x -= IFace.ivel.x;  Lft.vel.y -= IFace.ivel.y;  Lft.vel.z -= IFace.ivel.z;
    Rght.vel.x -= IFace.ivel.x; Rght.vel.y -= IFace.ivel.y; Rght.vel.z -= IFace.ivel.z;
    IFace.ivel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
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
    case FLUX_HLLC:
        hllc(Lft, Rght, IFace);
	break;
    default:
        throw std::runtime_error("Invalid flux calculator.");
    } // end switch

    // perform divergence cleaning of the magnetic field
    // implemented by modifying the flux of B according to the method of Dedner et al.
    if ( G.MHD ) {
	F.divB = 0.5*(Rght.B.x + Lft.B.x); 
	//NOTE: This is NOT the real divergence of B, used as an indicator according to the method of Dedner et al.
	if (G.div_clean) {
	    F.B.x += Lft.psi + 0.5*(Rght.psi - Lft.psi) - (G.c_h/2.0)*(Rght.B.x - Lft.B.x);
	    F.psi  = (F.divB - (1.0/(2.0*G.c_h))*(Rght.psi - Lft.psi)) * G.c_h * G.c_h;
	}
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
    // See Ian Johnston's thesis, Equations 3.31, 3.34 and 3.37
    /* Flux of Total Energy */
    F.total_energy += 0.5 * F.mass * pow(vabs(IFace.ivel),2) + dot(F.momentum, IFace.ivel);
    /* Flux of momentum */
    F.momentum += F.mass * IFace.ivel;
 
    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    // also transform the interface velocities
    IFace.ivel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
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
    fs.vel -= IFace.ivel; // Beware: fs.vel is changed here and restored below.
    IFace.ivel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
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
    F.total_energy += 0.5 * F.mass * pow(vabs(IFace.ivel),2) + dot(F.momentum, IFace.ivel);
    /* Flux of momentum */
    F.momentum += F.mass * IFace.ivel;

    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    // also transform the interface velocities
    IFace.ivel.transform_to_global(IFace.n, IFace.t1, IFace.t2);	  
    // also transform the magnetic field
    if ( G.MHD ) {
	F.B.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    }
    fs.vel.x = vx; fs.vel.y = vy; fs.vel.z = vz; // restore fs.vel
    return SUCCESS;
} // end compute_boundary_flux()

// added the artificial viscosity flux limiters
int artificial_diffusion(FV_Interface &IFace, FV_Cell &cL1, FV_Cell &cL0, FV_Cell &cR0, FV_Cell &cR1)
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();
    
    double rL1, rL0, rR0, rR1;
    double pL1, pL0, pR0, pR1;
    double aL0, aR0;    
    double uL1, uL0, uR0, uR1;
    double vL1, vL0, vR0, vR1;
    double wL1, wL0, wR0, wR1;
    double HL1, HL0, HR0, HR1;
    double pLrL1, pLrL0, pRrR0, pRrR1;
    double eL1, eL0, eR0, eR1;
    double keL1, keL0, keR0, keR1;
    
    // transfrom velocity and momentum fluxinto local frame of reference
    cL1.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    cL0.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    cR0.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    cR1.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
    
    IFace.ivel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
        
    IFace.F->momentum.transform_to_local(IFace.n, IFace.t1, IFace.t2);

    rL1 = cL1.fs->gas->rho;    rL0 = cL0.fs->gas->rho;
    pL1 = cL1.fs->gas->p;      pL0 = cL0.fs->gas->p;
    pLrL1 = pL1 / rL1;         pLrL0 = pL0 / rL0;
    uL1 = cL1.fs->vel.x;       uL0 = cL0.fs->vel.x;
    vL1 = cL1.fs->vel.y;       vL0 = cL0.fs->vel.y;
    wL1 = cL1.fs->vel.z;       wL0 = cL0.fs->vel.z;
    eL1 = accumulate(cL1.fs->gas->e.begin(), cL1.fs->gas->e.end(), 0.0);
    eL0 = accumulate(cL0.fs->gas->e.begin(), cL0.fs->gas->e.end(), 0.0);
    keL1 = 0.5 * (uL1 * uL1 + vL1 * vL1 + wL1 * wL1);
    keL0 = 0.5 * (uL0 * uL0 + vL0 * vL0 + wL0 * wL0);
    HL1 = eL1 + pLrL1 + keL1;      HL0 = eL0 + pLrL0 + keL0;   
    aL0 = cL0.fs->gas->a;        

    rR1 = cR1.fs->gas->rho;    rR0 = cR0.fs->gas->rho;
    pR1 = cR1.fs->gas->p;      pR0 = cR0.fs->gas->p;
    pRrR1 = pR1 / rR1;         pRrR0 = pR0 / rR0;
    uR1 = cR1.fs->vel.x;       uR0 = cR0.fs->vel.x;
    vR1 = cR1.fs->vel.y;       vR0 = cR0.fs->vel.y;
    wR1 = cR1.fs->vel.z;       wR0 = cR0.fs->vel.z;    
    eR1 = accumulate(cR1.fs->gas->e.begin(), cR1.fs->gas->e.end(), 0.0);
    eR0 = accumulate(cR0.fs->gas->e.begin(), cR0.fs->gas->e.end(), 0.0);
    keR1 = 0.5 * (uR1 * uR1 + vR1 * vR1 + wR1 * wR1);
    keR0 = 0.5 * (uR0 * uR0 + vR0 * vR0 + wR0 * wR0);
    HR1 = eR1 + pRrR1 + keR1;      HR0 = eR0 + pRrR0 + keR0;
    aR0 = cR0.fs->gas->a;                                       
    
    // variables
    double sensor_i, sensor_i1;
    double gamma_i, gamma_i1;    
    double epsilon_2, epsilon_4;
    double UL_1, UL_0, UR_0, UR_1;
    
    // constant
    double kappa_2 = G.artificial_kappa_2;
    double kappa_4 = G.artificial_kappa_4;
    
    // pressure sensor
    sensor_i = fabs( (pR0 - 2.0*pL0 + pL1) / (pR0 + 2.0*pL0 + pL1) );
    sensor_i1 = fabs( (pR1 - 2.0*pR0 + pL0) / (pR1 + 2.0*pR0 + pL0) );
    
    // spectral radius of the Jacobian matrix
    gamma_i = fabs( uL0-IFace.ivel.x ) + aL0;
    gamma_i1 = fabs( uR0-IFace.ivel.x ) + aR0;
    
    epsilon_2 = kappa_2 * max(sensor_i, sensor_i1) * max(gamma_i, gamma_i1);
    epsilon_4 = max(0.0, kappa_4*max(gamma_i, gamma_i1) - epsilon_2);  
    
    // mass
    UL_1 = rL1;
    UL_0 = rL0;
    UR_0 = rR0;
    UR_1 = rR1;
    IFace.F->mass -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );
        
    // momentum x
    UL_1 = rL1 * uL1;
    UL_0 = rL0 * uL0;
    UR_0 = rR0 * uR0;
    UR_1 = rR1 * uR1;
    IFace.F->momentum.x -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );
    
    // momentum y
    UL_1 = rL1 * vL1;
    UL_0 = rL0 * vL0;
    UR_0 = rR0 * vR0;
    UR_1 = rR1 * vR1;
    IFace.F->momentum.y -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );
    
    // momentum z
    UL_1 = rL1 * wL1;
    UL_0 = rL0 * wL0;
    UR_0 = rR0 * wR0;
    UR_1 = rR1 * wR1;        
    IFace.F->momentum.z -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );
    
    // total energy
    UL_1 = rL1 * HL1;
    UL_0 = rL0 * HL0;
    UR_0 = rR0 * HR0;
    UR_1 = rR1 * HR1;       
    IFace.F->total_energy -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );
    
    // tke
    UL_1 = rL1 * cL1.fs->tke;
    UL_0 = rL0 * cL0.fs->tke;
    UR_0 = rR0 * cR0.fs->tke;
    UR_1 = rR1 * cR1.fs->tke;       
    IFace.F->tke -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );
    
    // omega
    UL_1 = rL1 * cL1.fs->omega;
    UL_0 = rL0 * cL0.fs->omega;
    UR_0 = rR0 * cR0.fs->omega;
    UR_1 = rR1 * cR1.fs->omega;       
    IFace.F->omega -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );    
    
    // mass fraction
    for ( int isp = 0; isp < nsp; ++isp ) {
        UL_1 = rL1 * cL1.fs->gas->massf[isp];
        UL_0 = rL0 * cL0.fs->gas->massf[isp];
        UR_0 = rR0 * cR0.fs->gas->massf[isp];
        UR_1 = rR1 * cR1.fs->gas->massf[isp];
        IFace.F->massf[isp] -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );
    }
    
    // individual energies
    for ( int imode = 1; imode < nmodes; ++imode ) {
        UL_1 = rL1 * cL1.fs->gas->e[imode];
        UL_0 = rL0 * cL0.fs->gas->e[imode];
        UR_0 = rR0 * cR0.fs->gas->e[imode];
        UR_1 = rR1 * cR1.fs->gas->e[imode];    
	IFace.F->energies[imode] -= epsilon_2 * ( UR_0-UL_0 ) - epsilon_4 * ( UR_1-3.0*UR_0+3.0*UL_0-UL_1 );
    }   
    
    // transfrom velocity and momentum into global frame of reference
    cL1.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    cL0.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    cR0.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
    cR1.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2); 
    
    IFace.ivel.transform_to_global(IFace.n, IFace.t1, IFace.t2);    

    IFace.F->momentum.transform_to_global(IFace.n, IFace.t1, IFace.t2);        
   
    return SUCCESS;        
}
