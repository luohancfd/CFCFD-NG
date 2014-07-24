/**
 * fluxcalc.d
 * Flux calculators, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-23: initial cut, to explore options.
 */

module fluxcalc;

import std.math;
import std.stdio;
import std.conv;
import geom;
import gasmodelutil;
import flowstate;
import conservedquantities;
import fvcore;
import fvinterface;
import globalconfig;


void compute_interface_flux(ref FlowState Lft, ref FlowState Rght,
			    ref FVInterface IFace, double omegaz=0.0)
/** \brief Compute the inviscid fluxes (in 2D) across the cell interfaces.
 *
 * This is the top-level function that calls the previously selected
 * flux calculator.
 * Much of the detailed work is delegated.
 *
 * \param Lft    : reference to the LEFT flow state
 * \param Rght   : reference to the RIGHT flow state
 * \param IFace  : reference to the interface where the fluxes are to be stored
 *
 * Note that the FlowState objects are tampered with, but should be put back right
 * by the end of the function. [TODO] check this assertion.
 */
{
    double WSL, WSR;
    ConservedQuantities F = IFace.F;

    // Transform to interface frame of reference.
    Lft.vel -= IFace.vel;
    Rght.vel -= IFace.vel;
    IFace.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    Lft.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    Rght.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);

    // also transform the magnetic field
    if ( GlobalConfig.MHD ) {
	Lft.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        Rght.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    }

    // Compute the fluxes in the local frame of the interface.
    switch ( GlobalConfig.flux_calculator ) {
    case flux_efm:
        efmflx(Lft, Rght, IFace);
	break;
    case flux_ausmdv:
        ausmdv(Lft, Rght, IFace);
	break;
    case flux_adaptive:
        adaptive_flux(Lft, Rght, IFace);
	break;
    case flux_ausm_plus_up:
        ausm_plus_up(Lft, Rght, IFace);
	break;
    case flux_hlle:
        hlle(Lft, Rght, IFace);
	break;
    default:
        throw new Error(text("Invalid flux calculator, flux_claculator=.",
			     GlobalConfig.flux_calculator));
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
    
    // Transform fluxes back from interface frame of reference to local frame of reference.
    // Flux of Total Energy
    F.total_energy += 0.5 * F.mass * pow(abs(IFace.vel),2) + dot(F.momentum, IFace.vel);
    // Flux of momentum
    F.momentum += F.mass * IFace.vel;
 
    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    // also transform the interface velocities
    IFace.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    // also transform the magnetic field
    if ( GlobalConfig.MHD ) {
	F.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    }
	
    debug ( 2 ) {
	writeln("Interface fluxes");
	writeln("xyz_mom=", F.momentum);
	if ( GlobalConfig.MHD ) {
	    writeln("xyz_B=", F.B);
	}
    } // end debug
    
}

void set_flux_vector_in_local_frame(ref ConservedQuantities F, ref FlowState fs)
{
    double rho = fs.gas.rho;
    double un = fs.vel.x;
    double vt1 = fs.vel.y;
    double vt2 = fs.vel.z;
    double p = fs.gas.p;
    double e = 0.0; foreach(elem; fs.gas.e) e += elem;
    double ke = 0.5 * (un*un + vt1*vt1 + vt2*vt2); // Kinetic energy per unit volume.
    
    // Mass flux (mass / unit time / unit area)
    F.mass = rho * un; // The mass flux is relative to the moving interface.
    // Flux of normal momentum
    F.momentum.refx = F.mass * un + p;
    // Flux of tangential momentum
    F.momentum.refy = F.mass * vt1;
    F.momentum.refz = F.mass * vt2;
    // Flux of Total Energy
    F.total_energy = F.mass * (e + ke) + p * un;
    F.tke = F.mass * fs.tke;  // turbulence kinetic energy
    F.omega = F.mass * fs.omega;  // pseudo vorticity
    // Species mass flux
    for ( size_t isp = 0; isp < F.massf.length; ++isp ) {
	F.massf[isp] = F.mass * fs.gas.massf[isp];
    }
    // Individual energies.
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    for ( size_t imode = 1; imode < F.energies.length; ++imode ) {
	F.energies[imode] = F.mass * fs.gas.e[imode];
    }
}

void set_flux_vector_in_global_frame(ref FVInterface IFace, ref FlowState fs, 
				     double omegaz=0.0)
{
    ConservedQuantities F = IFace.F;
    // Record velocity to restore fs at end.
    double vx = fs.vel.x; double vy = fs.vel.y; double vz = fs.vel.z; 
    // Transform to interface frame of reference.
    fs.vel -= IFace.vel; // Beware: fs.vel is changed here and restored below.
    IFace.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    // also transform the magnetic field
    if ( GlobalConfig.MHD ) {
	fs.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    }
    set_flux_vector_in_local_frame(IFace.F, fs);
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
    F.total_energy += 0.5 * F.mass * pow(abs(IFace.vel),2) + dot(F.momentum, IFace.vel);
    /* Flux of momentum */
    F.momentum += F.mass * IFace.vel;

    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    // also transform the interface velocities
    IFace.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);	  
    // also transform the magnetic field
    if ( GlobalConfig.MHD ) {
	F.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    }
    fs.vel.refx = vx; fs.vel.refy = vy; fs.vel.refz = vz; // restore fs.vel
}

void ausmdv(in FlowState Lft, in FlowState Rght, ref FVInterface IFace)
{
}

void efmflx(in FlowState Lft, in FlowState Rght, ref FVInterface IFace)
/** \brief Compute the fluxes across an interface using
 * the Equilibrium Flux	Method of Macrossan & Pullin
 *
 * \param Lft    : IN     : array of Left flow states
 *     (with velocities in local frame of reference)
 * \param Rght   : IN     : array of Right flow state
 * \param IF     : IN/OUT : array of interface flux data structures
 *
 * \verbatim
 * interface data -- contains...
 *     flux of mass across the interface (kg/s/m**2)
 *     flux of normal momentum
 *     flux of tangential momentum
 *     flux of energy
 *     array of species fluxes 
 *     vibrational energies
 *     free-electron energy
 * \endverbatim 
 */
{
    /*
     * Local variable names reflect the names used in the original
     * FORTRAN code by MNM.
     */
    const double PHI = 1.0;
    double vnL, vpL, vnR, vpR, vqL, vqR;
    double rtL, cmpL, rtR, cmpR;
    double hvsqL, hvsqR;
    double wL, wR, dL, dR;
    double rhoL, rhoR, presL, presR, tR, tL;
    double eL, eR, hL, hR;
    double snL, snR, exL, exR, efL, efR;
    double fmsL, fmsR;
    double cv, cp, con, gam, Rgas;
    double cvL, cvR, RgasL, RgasR;
    double rLsqrt, rRsqrt, alpha;
    auto gmodel = GlobalConfig.gmodel;
    int statusf;

    /* Collect the global gas constants for later use. */
    int nsp = gmodel.n_species;
    int nmodes = gmodel.n_modes;

    /*
     * Calculate Constants
     */
    /* dtwspi = 1.0 / (2.0 * sqrt ( 3.14159265359 )); */
    const double dtwspi = 0.282094792;

    /*
     * Unpack Left flow state
     */
    rhoL = Lft.gas.rho;
    presL = Lft.gas.p;
    eL = 0.0; foreach(elem; Lft.gas.e) eL += elem;
    hL = eL + presL/rhoL;
    tL = Lft.gas.T[0];
    vnL = Lft.vel.x;
    vpL = Lft.vel.y;
    vqL = Lft.vel.z;

    /*
     * Unpack Right flow state
     */
    rhoR = Rght.gas.rho;
    presR = Rght.gas.p;
    eR = 0.0; foreach(elem; Rght.gas.e) eR += elem;
    hR = eR + presR/rhoR;
    tR = Rght.gas.T[0];
    vnR = Rght.vel.x;
    vpR = Rght.vel.y;
    vqR = Rght.vel.z;

    /* Derive the gas "constants" from the local conditions. */
    cvL = gmodel.Cv(Lft.gas);
    RgasL = presL / (rhoL * tL);
    cvR = gmodel.Cv(Rght.gas);
    RgasR = presR / (rhoR * tR);

    rLsqrt = sqrt(rhoL);
    rRsqrt = sqrt(rhoR);
    alpha = rLsqrt / (rLsqrt + rRsqrt);

    cv = alpha * cvL + (1.0 - alpha) * cvR;
    Rgas = alpha * RgasL + (1.0 - alpha) * RgasR;

    cp = cv + Rgas;
    gam = cp / cv;

    /*
     * Start EFM calculation proper.
     */
    con = 0.5 * (gam + 1.0) / (gam - 1.0);

    rtL = Rgas * tL;
    cmpL = sqrt(2.0 * rtL);
    hvsqL = 0.5 * (vnL * vnL + vpL * vpL + vqL * vqL);

    snL = vnL / (PHI * cmpL);
    exxef(snL, exL, efL);

    wL = 0.5 * (1.0 + efL);
    dL = exL * dtwspi;

    rtR = presR / rhoR;
    cmpR = sqrt(2.0 * rtR);
    hvsqR = 0.5 * (vnR * vnR + vpR * vpR + vqR * vqR);

    snR = vnR / (PHI * cmpR);
    exxef(snR, exR, efR);

    wR = 0.5 * (1.0 - efR);
    dR = -exR * dtwspi;

    /*
     * combine the fluxes
     */
    fmsL = (wL * rhoL * vnL) + (dL * cmpL * rhoL);
    fmsR = (wR * rhoR * vnR) + (dR * cmpR * rhoR);

    ConservedQuantities F = IFace.F;
    F.mass = fmsL + fmsR;

    F.momentum.refx = fmsL * vnL + fmsR * vnR + wL * presL + wR * presR;
    F.momentum.refy = fmsL * vpL + fmsR * vpR;
    F.momentum.refz = fmsL * vqL + fmsR * vqR;

    F.total_energy = (wL * rhoL * vnL) * (hvsqL + hL) +
	(wR * rhoR * vnR) * (hvsqR + hR) +
	(dL * cmpL * rhoL) * (hvsqL + con * rtL) +
	(dR * cmpR * rhoR) * (hvsqR + con * rtR);

    if (F.mass > 0.0) {
	F.tke = F.mass * Lft.tke;
	F.omega = F.mass * Lft.omega;
    } else {
	F.tke = F.mass * Rght.tke;
	F.omega = F.mass * Rght.omega;
    }

    /* 
     * Species mass flux.
     * Presently, this is implemented by assuming that
     * the wind is blowing one way or the other and then
     * picking the appropriate side for the species fractions.
     * Such an approach may not be fully compatible with the
     * EFM approach where there can be fluxes from both sides.
     */
    for ( int isp = 0; isp < nsp; ++isp ) {
	if (F.mass > 0.0)
	    F.massf[isp] = (F.mass) * Lft.gas.massf[isp];
	else
	    F.massf[isp] = (F.mass) * Rght.gas.massf[isp];
    }   /* isp loop */

    // Individual energies.
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    if (F.mass > 0.0) {
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = (F.mass) * Lft.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += (F.mass) * Lft.gas.p_e / Lft.gas.rho;
    } else {
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = F.mass * Rght.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += F.mass * Rght.gas.p_e / Rght.gas.rho;
    }
} // end efmflx()

void exxef(double sn, ref double exx, ref double ef)
/** \brief Compute exp(-x**2) and erf(x) with a polynomial approximation.
 *
 * \param sn   : IN  : x
 * \param &exx : OUT : exp(x**2)
 * \param &ef  : OUT : erf(x)  error function
 */
{
    double snsq, ef1, y;

    const double P = 0.327591100;
    const double A1 = 0.254829592;
    const double A2 = -0.284496736;
    const double A3 = 1.421413741;
    const double A4 = -1.453152027;
    const double A5 = 1.061405429;
    const double LIMIT = 5.0;
    const double EXLIM = 0.138879e-10;
    const double EFLIM = 1.0;

    //#   define DSIGN(val,sgn) ( (sgn >= 0.0)? fabs(val): -fabs(val) )

    if (fabs(sn) > LIMIT) {
        exx = EXLIM;
        ef1 = EFLIM;
    } else {
        snsq = sn * sn;
        exx = exp(-snsq);
        y = 1.0 / (1.0 + P * fabs(sn));
        ef1 = 1.0 - y * (A1 + y * (A2 + y * (A3 + y * (A4 + A5 * y)))) * exx;
    }
    ef = copysign(ef1, sn);
} // end exxef

void adaptive_flux(in FlowState Lft, in FlowState Rght, ref FVInterface IFace)
{
}

void ausm_plus_up(in FlowState Lft, in FlowState Rght, ref FVInterface IFace)
{
}

void hlle(in FlowState Lft, in FlowState Rght, ref FVInterface IFace)
{
}
