/// \file one_d_interp.cxx

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
// #include <iostream>
// #include <string>
// #include <sstream>

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/geometry2/source/geom.hh"
#include "cell.hh"
#include "one_d_interp.hh"
#include "kernel.hh"

thermo_interp_t thermo_interpolator = INTERP_RHOE;

thermo_interp_t set_thermo_interpolator(thermo_interp_t interp)
{
    thermo_interpolator = interp;
    return thermo_interpolator;
}

thermo_interp_t get_thermo_interpolator()
{
    return thermo_interpolator;
}

std::string get_thermo_interpolator_name(thermo_interp_t interp)
{
    switch ( interp ) {
    case INTERP_PT: return "pT";
    case INTERP_RHOE: return "rhoe";
    case INTERP_RHOP: return "rhop";
    case INTERP_RHOT: return "rhoT";
    default: return "none";
    }
} // end get_thermo_interpolator_name()

//-----------------------------------------------------------------------------

/// \brief One-dimensional reconstruction of a scalar quantity.
///
/// See MBCNS workbook 2000/2 page 36 (26-Jan-2001) for formulation.
/// and MBCNS workbook 2005/Apr page 36 for new index labels
///
inline int one_d_interp_scalar( double qL1, double qL0, double qR0, double qR1, 
				double lenL1, double lenL0, double lenR0, double lenR1, 
				double &qL, double &qR, int apply_limiter,
				int extrema_clipping )
{
    double aL0, aR0, delLminus, del, delRplus, sL, sR;
    double lower_limit, upper_limit;
    // Set up differences and limiter values.
    aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
    aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
    delLminus = 2.0 * (qL0 - qL1) / (lenL0 + lenL1);
    del = 2.0 * (qR0 - qL0) / (lenR0 + lenL0);
    delRplus = 2.0 * (qR1 - qR0) / (lenR1 + lenR0);
    if ( apply_limiter == 1 ) {
	// val Albada limiter as per Ian Johnston's thesis.
	const double epsilon = 1.0e-12;
	sL = (delLminus * del + fabs(delLminus * del)) /
	    (delLminus * delLminus + del * del + epsilon);
	sR = (del * delRplus + fabs(del * delRplus)) /
	    (del * del + delRplus * delRplus + epsilon);
    } else {
	sL = 1.0;
	sR = 1.0;
    }
    // The high-order reconstruction, possibly limited.
    qL = qL0 + sL * aL0 * ( del * (2.0*lenL0 + lenL1) + delLminus * lenR0 );
    qR = qR0 - sR * aR0 * ( delRplus * lenL0 + del * (2.0*lenR0 + lenR1) );
    if ( extrema_clipping == 1 ) {
	// An extra limiting filter to make sure that we have not introduced
	// any new extreme values.
	// This was introduced to deal with very sharp transitions in species.
	lower_limit = MINIMUM(qL0, qR0);
	upper_limit = MAXIMUM(qL0, qR0);
	qL = MINIMUM(upper_limit, MAXIMUM(lower_limit, qL)); 
	qR = MINIMUM(upper_limit, MAXIMUM(lower_limit, qR));
    }
    return SUCCESS;
} // end of one_d_interp_scalar()

//-----------------------------------------------------------------------------

/// \brief Reconstruct flow properties at an interface from FV_Cell properties.
///
/// This is essentially a one-dimensional interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.
int one_d_interp(const FV_Cell &cL1, const FV_Cell &cL0,
		 const FV_Cell &cR0, const FV_Cell &cR1,
		 double cL1Length, double cL0Length,
		 double cR0Length, double cR1Length,
		 FlowState &Lft, FlowState &Rght )
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();

    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Lft.copy_values_from(*(cL0.fs));
    Rght.copy_values_from(*(cR0.fs));
    if ( get_Xorder_flag() > 1 ) {
	// High-order reconstruction for some properties.
	one_d_interp_scalar(cL1.fs->vel.x, cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x,
			    cL1Length, cL0Length, cR0Length, cR1Length,
			    Lft.vel.x, Rght.vel.x, apply_limiter_flag, extrema_clipping_flag);
	one_d_interp_scalar(cL1.fs->vel.y, cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y,
			    cL1Length, cL0Length, cR0Length, cR1Length,
			    Lft.vel.y, Rght.vel.y, apply_limiter_flag, extrema_clipping_flag);
	one_d_interp_scalar(cL1.fs->vel.z, cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z,
			    cL1Length, cL0Length, cR0Length, cR1Length,
			    Lft.vel.z, Rght.vel.z, apply_limiter_flag, extrema_clipping_flag);
	if ( get_mhd_flag() == 1 ) {
	    one_d_interp_scalar(cL1.fs->B.x, cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.B.x, Rght.B.x, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interp_scalar(cL1.fs->B.y, cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.B.y, Rght.B.y, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interp_scalar(cL1.fs->B.z, cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.B.z, Rght.B.z, apply_limiter_flag, extrema_clipping_flag);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    one_d_interp_scalar(cL1.fs->tke, cL0.fs->tke, cR0.fs->tke, cR1.fs->tke,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.tke, Rght.tke, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interp_scalar(cL1.fs->omega, cL0.fs->omega, cR0.fs->omega, cR1.fs->omega,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.omega, Rght.omega, apply_limiter_flag, extrema_clipping_flag);
	}
	Gas_data &gL1 = *(cL1.fs->gas);
	Gas_data &gL0 = *(cL0.fs->gas);
	Gas_data &gR0 = *(cR0.fs->gas);
	Gas_data &gR1 = *(cR1.fs->gas);
	if ( nsp > 1 ) {
	    // Multiple species.
	    for ( size_t isp = 0; isp < nsp; ++isp ) {
		one_d_interp_scalar(gL1.massf[isp], gL0.massf[isp], gR0.massf[isp], gR1.massf[isp],
				    cL1Length, cL0Length, cR0Length, cR1Length,
				    Lft.gas->massf[isp], Rght.gas->massf[isp],
				    apply_limiter_flag, extrema_clipping_flag);
	    }
	    if ( scale_mass_fractions( Lft.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp < nsp; ++isp )
		    Lft.gas->massf[isp] = cL0.fs->gas->massf[isp];
	    }
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp < nsp; ++isp )
		    Rght.gas->massf[isp] = cR0.fs->gas->massf[isp];
	    }
	} else {
	    // Only one possible mass-fraction value for a single species.
	    Lft.gas->massf[0] = 1.0;
	    Rght.gas->massf[0] = 1.0;
	}

	// Interpolate on two of the thermodynamic quantities, 
	// and fill in the rest based on an EOS call. 
	// If an EOS call fails, fall back to just copying cell-centre data.
	// This does presume that the cell-centre data is valid. 
	switch ( thermo_interpolator ) {
	case INTERP_PT: 
	    one_d_interp_scalar(gL1.p, gL0.p, gR0.p, gR1.p,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.gas->p, Rght.gas->p, apply_limiter_flag, extrema_clipping_flag);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_scalar(gL1.T[i], gL0.T[i], gR0.T[i], gR1.T[i],
				    cL1Length, cL0Length, cR0Length, cR1Length,
				    Lft.gas->T[i], Rght.gas->T[i], apply_limiter_flag, extrema_clipping_flag);
	    }
	    if ( gmodel->eval_thermo_state_pT(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    if ( gmodel->eval_thermo_state_pT(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOE:
	    one_d_interp_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.gas->rho, Rght.gas->rho, apply_limiter_flag, extrema_clipping_flag);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_scalar(gL1.e[i], gL0.e[i], gR0.e[i], gR1.e[i],
				    cL1Length, cL0Length, cR0Length, cR1Length,
				    Lft.gas->e[i], Rght.gas->e[i], apply_limiter_flag, extrema_clipping_flag);
	    }
	    if ( gmodel->eval_thermo_state_rhoe(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    if ( gmodel->eval_thermo_state_rhoe(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOP:
	    one_d_interp_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.gas->rho, Rght.gas->rho, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interp_scalar(gL1.p, gL0.p, gR0.p, gR1.p,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.gas->p, Rght.gas->p, apply_limiter_flag, extrema_clipping_flag);
	    if ( gmodel->eval_thermo_state_rhop(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    if ( gmodel->eval_thermo_state_rhop(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOT: 
	    one_d_interp_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.gas->rho, Rght.gas->rho, apply_limiter_flag, extrema_clipping_flag);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_scalar(gL1.T[i], gL0.T[i], gR0.T[i], gR1.T[i],
				    cL1Length, cL0Length, cR0Length, cR1Length,
				    Lft.gas->T[i], Rght.gas->T[i], apply_limiter_flag, extrema_clipping_flag);
	    }
	    if ( gmodel->eval_thermo_state_rhoT(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    if ( gmodel->eval_thermo_state_rhoT(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	default: 
	    throw std::runtime_error("Invalid thermo interpolator.");
	}

	if ( get_viscous_flag() ) {
	    // FIX-ME -- for speed, just copy the cell data.
	    gmodel->eval_transport_coefficients(*(Lft.gas));
	    gmodel->eval_transport_coefficients(*(Rght.gas));
	}
	if ( get_diffusion_flag() ) {
	    // FIX-ME -- for speed, just copy the cell data.
	    gmodel->eval_diffusion_coefficients(*(Lft.gas));
	    gmodel->eval_diffusion_coefficients(*(Rght.gas));
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of one_d_interp()

//-----------------------------------------------------------------------------

inline int mach_weighted_one_d_interp_scalar( double qL1, double qL0, double qR0, double qR1, 
				 double lenL1, double lenL0, double lenR0, double lenR1, 
				 double &qL, double &qR, double kL, double kR,
				 int apply_limiter, int extrema_clipping )
{
    double aL0, aR0, delLminus, del, delRplus, sL, sR;
    double lower_limit, upper_limit;
    // Set up differences and limiter values.
    aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
    aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
    delLminus = 2.0 * (qL0 - qL1) / (lenL0 + lenL1);
    del = 2.0 * (qR0 - qL0) / (lenR0 + lenL0);
    delRplus = 2.0 * (qR1 - qR0) / (lenR1 + lenR0);
    if ( apply_limiter == 1 ) {
	// val Albada limiter as per Ian Johnston's thesis.
	const double epsilon = 1.0e-12;
	sL = (delLminus * del + fabs(delLminus * del)) /
	    (delLminus * delLminus + del * del + epsilon);
	sR = (del * delRplus + fabs(del * delRplus)) /
	    (del * del + delRplus * delRplus + epsilon);
    } else {
	sL = 1.0;
	sR = 1.0;
    }
    
    // The high-order reconstruction, possibly limited.
    qL = qL0 + sL * aL0 * ( (1 + kL*sL) * del * (2.0*lenL0 + lenL1) + (1 - kL*sL) * delLminus * lenR0 );
    qR = qR0 - sR * aR0 * ( (1 - kR*sR) * delRplus * lenL0 + (1 + kR*sR) * del * (2.0*lenR0 + lenR1) );
    if ( extrema_clipping == 1 ) {
	// An extra limiting filter to make sure that we have not introduced
	// any new extreme values.
	// This was introduced to deal with very sharp transitions in species.
	lower_limit = MINIMUM(qL0, qR0);
	upper_limit = MAXIMUM(qL0, qR0);
	qL = MINIMUM(upper_limit, MAXIMUM(lower_limit, qL)); 
	qR = MINIMUM(upper_limit, MAXIMUM(lower_limit, qR));
    }
    return SUCCESS;
} // end of mach_weighted_one_d_interp_scalar()


/// \brief Reconstruct upstream Mach number weighted flow properties at an interface from FV_Cell properties.
///
/// This is essentially a one-dimensional interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.
int mach_weighted_one_d_interp(const FV_Cell &cL1, const FV_Cell &cL0,
			       const FV_Cell &cR0, const FV_Cell &cR1,
			       double cL1Length, double cL0Length,
			       double cR0Length, double cR1Length,
			       FlowState &Lft, FlowState &Rght )
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data &gL1 = *(cL1.fs->gas);
    Gas_data &gL0 = *(cL0.fs->gas);
    Gas_data &gR0 = *(cR0.fs->gas);
    Gas_data &gR1 = *(cR1.fs->gas);
    int nsp = gmodel->get_number_of_species();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();

    // If flow is supersonic in the direction of the interface don't use any
    // downstream information. Else, linearly velocity weight the upstream information.
    double ML = dot(cL0.fs->vel, unit(cR0.pos[0] - cL0.pos[0])) / cL0.fs->gas->a;
    double MR = dot(cR0.fs->vel, unit(cL0.pos[0] - cR0.pos[0])) / cR0.fs->gas->a;
    double kL = min(1.0, max(0.0, ML));
    double kR = min(1.0, max(0.0, MR));

    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Lft.copy_values_from(*(cL0.fs));
    Rght.copy_values_from(*(cR0.fs));
    if ( get_Xorder_flag() > 1 ) {
	// High-order reconstruction for some properties.
	mach_weighted_one_d_interp_scalar(cL1.fs->vel.x, cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x,
					  cL1Length, cL0Length, cR0Length, cR1Length,
					  Lft.vel.x, Rght.vel.x, kL, kR, apply_limiter_flag,
					  extrema_clipping_flag);
	mach_weighted_one_d_interp_scalar(cL1.fs->vel.y, cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y,
					  cL1Length, cL0Length, cR0Length, cR1Length,
					  Lft.vel.y, Rght.vel.y, kL, kR, apply_limiter_flag,
					  extrema_clipping_flag);
	mach_weighted_one_d_interp_scalar(cL1.fs->vel.z, cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z,
					  cL1Length, cL0Length, cR0Length, cR1Length,
					  Lft.vel.z, Rght.vel.z, kL, kR, apply_limiter_flag, 
					  extrema_clipping_flag);
	if ( get_mhd_flag() == 1 ) {
	    mach_weighted_one_d_interp_scalar(cL1.fs->B.x, cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x,
					      cL1Length, cL0Length, cR0Length, cR1Length,
					      Lft.B.x, Rght.B.x, kL, kR, apply_limiter_flag, 
					      extrema_clipping_flag);
	    mach_weighted_one_d_interp_scalar(cL1.fs->B.y, cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y,
					      cL1Length, cL0Length, cR0Length, cR1Length,
					      Lft.B.y, Rght.B.y, kL, kR, apply_limiter_flag,
					      extrema_clipping_flag);
	    mach_weighted_one_d_interp_scalar(cL1.fs->B.z, cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z,
					      cL1Length, cL0Length, cR0Length, cR1Length,
					      Lft.B.z, Rght.B.z, kL, kR, apply_limiter_flag,
					      extrema_clipping_flag);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    mach_weighted_one_d_interp_scalar(cL1.fs->tke, cL0.fs->tke, cR0.fs->tke, cR1.fs->tke,
					      cL1Length, cL0Length, cR0Length, cR1Length,
					      Lft.tke, Rght.tke, kL, kR, apply_limiter_flag,
					      extrema_clipping_flag);
	    mach_weighted_one_d_interp_scalar(cL1.fs->omega, cL0.fs->omega, cR0.fs->omega, cR1.fs->omega,
					      cL1Length, cL0Length, cR0Length, cR1Length,
					      Lft.omega, Rght.omega, kL, kR, apply_limiter_flag, 
					      extrema_clipping_flag);
	}
	for ( int isp = 0; isp < nsp; ++isp ) {
	    mach_weighted_one_d_interp_scalar(cL1.fs->gas->massf[isp], cL0.fs->gas->massf[isp],
					      cR0.fs->gas->massf[isp], cR1.fs->gas->massf[isp],
					      cL1Length, cL0Length, cR0Length, cR1Length,
					      Lft.gas->massf[isp], Rght.gas->massf[isp],  kL, kR,
					      apply_limiter_flag, extrema_clipping_flag);
	}
	
	// Make the thermodynamic properties consistent.
	// Pressure, Local Speed of Sound and Temperature.
	// The value of 1 indicates that the old temperature
	// should be used as an initial guess for the iterative
	// EOS functions.
	// If the EOS call fouls up, just copy the cell data, low-order.
	if ( nsp > 1 ) {
	    if ( scale_mass_fractions( Lft.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Lft.gas->massf.size(); ++isp )
		    cL0.fs->gas->massf[isp] = Lft.gas->massf[isp];
	    }
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Rght.gas->massf.size(); ++isp )
		    cR0.fs->gas->massf[isp] = Rght.gas->massf[isp];
	    }
	}
	// Interpolate on two of the thermodynamic quantities, and fill
	// in the rest based on an EOS call.
	
	mach_weighted_one_d_interp_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho,
					  cL1Length, cL0Length, cR0Length, cR1Length,
					  Lft.gas->rho, Rght.gas->rho,
					  kL, kR,
					  apply_limiter_flag, extrema_clipping_flag);
	for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
	    mach_weighted_one_d_interp_scalar(gL1.e[i], gL0.e[i], gR0.e[i], gR1.e[i],
					      cL1Length, cL0Length, cR0Length, cR1Length,
					      Lft.gas->e[i], Rght.gas->e[i],
					      kL, kR,
					      apply_limiter_flag, extrema_clipping_flag);
	}

	if ( gmodel->eval_thermo_state_rhoe(*(Lft.gas)) != SUCCESS ) {
	    // Lft state failed, just copy values from adjacent cell
	    Lft.copy_values_from(*(cL0.fs));
	}
	if ( gmodel->eval_thermo_state_rhoe(*(Rght.gas)) != SUCCESS ) {
	    // Rght state failed, just copy values from adjacent cell
	    Rght.copy_values_from(*(cR0.fs));
	}
	
	if ( get_viscous_flag() ) {
	    gmodel->eval_transport_coefficients(*(Lft.gas));
	    gmodel->eval_transport_coefficients(*(Rght.gas));
	}
	if ( get_diffusion_flag() ) {
	    gmodel->eval_diffusion_coefficients(*(Lft.gas));
	    gmodel->eval_diffusion_coefficients(*(Rght.gas));
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of mach_weighted_one_d_interp()

//---------------------------------------------------------------------------


/// \brief One-sided one-dimensional reconstruction of a scalar quantity.
///
/// See Ian Johnston's thesis.
///
inline int onesided_interp_scalar( double qR0, double qR1, double qR2, 
				   double lenR0, double lenR1, double lenR2, 
				   double &qR, int apply_limiter,
				   int extrema_clipping )
{
    const double epsilon = 1.0e-12;

    double aR0, delRminus, delRplus, sR;
    
    // Set up differences and limiter values.
    aR0 = 0.5 * lenR0 / (lenR0 + 2.0*lenR1 + lenR2);
    delRminus = 2.0 * (qR1 - qR0) / (lenR1 + lenR0);
    delRplus = 2.0 * (qR2 - qR1) / (lenR2 + lenR1);
    if ( apply_limiter == 1 ) {
	// val Albada limiter as per Ian Johnston's thesis.
	sR = (delRminus * delRplus + fabs(delRminus * delRplus)) /
	    (delRminus * delRminus + delRplus * delRplus + epsilon);
    } else {
	sR = 1.0;
    }
    // The high-order reconstruction, possibly limited.
    qR = qR0 - 0.5 * aR0 * sR * ( (1 - sR) * delRplus * lenR1 + (1 + sR) * 
				  delRminus * (lenR1 + lenR2 + lenR0) );

    return SUCCESS;
} // end of onesided_interp_scalar()

/// \brief Asymmetric one-dimensional reconstruction of a scalar quantity.
///
/// See Ian Johnston's thesis.
///
inline int one_d_interior_interp_scalar( double qL0, double qR0, double qR1, double qR2, 
					 double lenL0, double lenR0, double lenR1, double lenR2,
					 double &qL, double &qR, 
					 int apply_limiter, int extrema_clipping )
{
    const double epsilon = 1.0e-12;

    double aL0, aR0, del, delRplus, delRminus, sL, sR;
    double lower_limit, upper_limit;
    
    // Set up differences and limiter values.
    aL0 = 0.5 * lenL0 / (lenL0 + 2.0*lenR0 + lenR1);
    aR0 = 0.5 * lenR0 / (lenR0 + 2.0*lenR1 + lenR2);
    del = 2.0 * (qR0 - qL0) / (lenR0 + lenL0);
    delRminus = 2.0 * (qR1 - qR0) / (lenR1 + lenR0);
    delRplus = 2.0 * (qR2 - qR1) / (lenR2 + lenR1);
    if ( apply_limiter == 1 ) {
	// val Albada limiter as per Ian Johnston's thesis.
	sL = (del * delRminus + fabs(del * delRminus)) /
	    (del * del + delRminus * delRminus + epsilon);
	sR = (delRminus * delRplus + fabs(delRminus * delRplus)) /
	    (delRminus * delRminus + delRplus * delRplus + epsilon);
    } else {
	sL = 1.0;
	sR = 1.0;
    }
    // The high-order reconstruction, possibly limited.
    qL = qL0 + aL0 * sL * ( delRminus * (2 * lenL0 + lenR0) + del * (lenR0 + lenR1 - lenL0) );
    qR = qR0 - aR0 * sR * ( delRplus * lenR1 + delRminus * (lenR1 + lenR2 + lenR0) );
    if ( extrema_clipping == 1 ) {
	// An extra limiting filter to make sure that we have not introduced
	// any new extreme values.
	// This was introduced to deal with very sharp transitions in species.
	lower_limit = MINIMUM(qL0, qR0);
	upper_limit = MAXIMUM(qL0, qR0);
	qL = MINIMUM(upper_limit, MAXIMUM(lower_limit, qL));
	qR = MINIMUM(upper_limit, MAXIMUM(lower_limit, qR));
    }
    return SUCCESS;
} // end of one_d_interior_interp_scalar()


/// \brief Reconstruct flow properties at an interface from FV_Cell properties from one side only.
///
/// This scheme uses the three cells to the right of the interface to do
/// a one-sided extrapolation of the flow properties. This is done to determine the shock
/// speed when shock-fitting.
int onesided_interp(const FV_Cell &cL0, const FV_Cell &cR0, const FV_Cell &cR1,
		    double cL0Length, double cR0Length, double cR1Length,
		    FlowState &Rght )
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data &gL0 = *(cL0.fs->gas);
    Gas_data &gR0 = *(cR0.fs->gas);
    Gas_data &gR1 = *(cR1.fs->gas);
    int nsp = gmodel->get_number_of_species();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();
    
    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Rght.copy_values_from(*(cR0.fs));
    if ( get_Xorder_flag() > 1 ) {
	// High-order reconstruction for some properties.
	onesided_interp_scalar(cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x,
			       cL0Length, cR0Length, cR1Length,
			       Rght.vel.x, apply_limiter_flag, extrema_clipping_flag);
	onesided_interp_scalar(cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y,
			       cL0Length, cR0Length, cR1Length,
			       Rght.vel.y, apply_limiter_flag, extrema_clipping_flag);
	onesided_interp_scalar(cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z,
			       cL0Length, cR0Length, cR1Length,
			       Rght.vel.z, apply_limiter_flag, extrema_clipping_flag);
	if ( get_mhd_flag() == 1 ) {
	    onesided_interp_scalar(cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x,
				   cL0Length, cR0Length, cR1Length,
				   Rght.B.x, apply_limiter_flag, extrema_clipping_flag);
	    onesided_interp_scalar(cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y,
				   cL0Length, cR0Length, cR1Length,
				   Rght.B.y, apply_limiter_flag, extrema_clipping_flag);
	    onesided_interp_scalar(cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z,
				   cL0Length, cR0Length, cR1Length,
				   Rght.B.z, apply_limiter_flag, extrema_clipping_flag);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    onesided_interp_scalar(cL0.fs->tke, cR0.fs->tke, cR1.fs->tke,
				   cL0Length, cR0Length, cR1Length,
				   Rght.tke, apply_limiter_flag, extrema_clipping_flag);
	    onesided_interp_scalar(cL0.fs->omega, cR0.fs->omega, cR1.fs->omega,
				   cL0Length, cR0Length, cR1Length,
				   Rght.omega, apply_limiter_flag, extrema_clipping_flag);
	}
        for ( int isp = 0; isp < nsp; ++isp ) {
	    onesided_interp_scalar(cL0.fs->gas->massf[isp],
				   cR0.fs->gas->massf[isp], cR1.fs->gas->massf[isp],
				   cL0Length, cR0Length, cR1Length,
				   Rght.gas->massf[isp],
				   apply_limiter_flag, extrema_clipping_flag);
        }
	
	// Make the thermodynamic properties consistent.
	// Pressure, Local Speed of Sound and Temperature.
        // The value of 1 indicates that the old temperature
        // should be used as an initial guess for the iterative
        // EOS functions.
	// If the EOS call fouls up, just copy the cell data, low-order.
	if ( nsp > 1 ) {
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Rght.gas->massf.size(); ++isp )
		    cR0.fs->gas->massf[isp] = Rght.gas->massf[isp];
	    }
	}
	
	// Interpolate on two of the thermodynamic quantities, and fill
	// in the rest based on an EOS call.
	
	onesided_interp_scalar(gL0.rho, gR0.rho, gR1.rho,
	                       cL0Length, cR0Length, cR1Length,
	                       Rght.gas->rho, apply_limiter_flag, extrema_clipping_flag);
	for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
	    onesided_interp_scalar(gL0.e[i], gR0.e[i], gR1.e[i],
				   cL0Length, cR0Length, cR1Length,
				   Rght.gas->e[i], apply_limiter_flag, extrema_clipping_flag);
	}
	int status = gmodel->eval_thermo_state_rhoe(*(Rght.gas));
	
	if ( status != SUCCESS ) {
	    // Rght state failed.
	    Rght.copy_values_from(*(cL0.fs));
	}			      
	
	if ( get_viscous_flag() ) {
	    gmodel->eval_transport_coefficients(*(Rght.gas));
	}
	if ( get_diffusion_flag() ) {
	    gmodel->eval_diffusion_coefficients(*(Rght.gas));
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of onesided_interp()

/// \brief Reconstruct flow properties at an interface from FV_Cell properties.
///
/// This scheme uses the three cells to the right of the interface to do
/// a one-sided extrapolation of the flow properties and two cells to the left and one to the right
/// to perform an interpolation of the flow properties. This is to maintain the same interpolation order 
/// when reconstructing the flow properties of the first interface after the shock when shock-fitting.
int one_d_interior_interp(const FV_Cell &cL0, const FV_Cell &cR0, 
			  const FV_Cell &cR1, const FV_Cell &cR2,
			  double cL0Length, double cR0Length,
			  double cR1Length, double cR2Length,
			  FlowState &Lft, FlowState &Rght )
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data &gL0 = *(cL0.fs->gas);
    Gas_data &gR0 = *(cR0.fs->gas);
    Gas_data &gR1 = *(cR1.fs->gas);
    Gas_data &gR2 = *(cR2.fs->gas);
    int nsp = gmodel->get_number_of_species();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();

    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Lft.copy_values_from(*(cL0.fs));
    Rght.copy_values_from(*(cR0.fs));
    if ( get_Xorder_flag() > 1 ) {
	// High-order reconstruction for some properties.
	one_d_interior_interp_scalar(cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x, cR2.fs->vel.x,
				     cL0Length, cR0Length, cR1Length, cR2Length,
				     Lft.vel.x, Rght.vel.x, apply_limiter_flag, extrema_clipping_flag);
	one_d_interior_interp_scalar(cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y, cR2.fs->vel.y,
				     cL0Length, cR0Length, cR1Length, cR2Length,
				     Lft.vel.y, Rght.vel.y, apply_limiter_flag, extrema_clipping_flag);
	one_d_interior_interp_scalar(cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z, cR2.fs->vel.z,
				     cL0Length, cR0Length, cR1Length, cR2Length,
				     Lft.vel.z, Rght.vel.z, apply_limiter_flag, extrema_clipping_flag);
	if ( get_mhd_flag() == 1 ) {
	    one_d_interior_interp_scalar(cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x, cR2.fs->B.x,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.B.x, Rght.B.x, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interior_interp_scalar(cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y, cR2.fs->B.y,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.B.y, Rght.B.y, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interior_interp_scalar(cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z, cR2.fs->B.z,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.B.z, Rght.B.z, apply_limiter_flag, extrema_clipping_flag);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    one_d_interior_interp_scalar(cL0.fs->tke, cR0.fs->tke, cR1.fs->tke, cR2.fs->tke,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.tke, Rght.tke, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interior_interp_scalar(cL0.fs->omega, cR0.fs->omega, cR1.fs->omega, cR2.fs->omega,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.omega, Rght.omega, apply_limiter_flag, extrema_clipping_flag);
	}
        for ( int isp = 0; isp < nsp; ++isp ) {
	    one_d_interior_interp_scalar(cL0.fs->gas->massf[isp], cR0.fs->gas->massf[isp],
					 cR1.fs->gas->massf[isp], cR2.fs->gas->massf[isp],
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.gas->massf[isp], Rght.gas->massf[isp],
					 apply_limiter_flag, extrema_clipping_flag);
        }

	// Make the thermodynamic properties consistent.
	// Pressure, Local Speed of Sound and Temperature.
        // The value of 1 indicates that the old temperature
        // should be used as an initial guess for the iterative
        // EOS functions.
	// If the EOS call fouls up, just copy the cell data, low-order.
	if ( nsp > 1 ) {
	    if ( scale_mass_fractions( Lft.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Lft.gas->massf.size(); ++isp )
		    cL0.fs->gas->massf[isp] = Lft.gas->massf[isp];
	    }
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Rght.gas->massf.size(); ++isp )
		    cR0.fs->gas->massf[isp] = Rght.gas->massf[isp];
	    }
	}

	// Interpolate on two of the thermodynamic quantities, and fill
	// in the rest based on an EOS call.
	one_d_interior_interp_scalar(gL0.rho, gR0.rho, gR1.rho, gR2.rho,
				     cL0Length, cR0Length, cR1Length, cR2Length,
				     Lft.gas->rho, Rght.gas->rho,
				     apply_limiter_flag, extrema_clipping_flag);

	for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
	    one_d_interior_interp_scalar(gL0.e[i], gR0.e[i], gR1.e[i], gR2.e[i],
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.gas->e[i], Rght.gas->e[i],
					 apply_limiter_flag, extrema_clipping_flag);
	}

	if ( get_viscous_flag() ) {
	    gmodel->eval_transport_coefficients(*(Lft.gas));
	    gmodel->eval_transport_coefficients(*(Rght.gas));
	}
	if ( get_diffusion_flag() ) {
	    gmodel->eval_diffusion_coefficients(*(Lft.gas));
	    gmodel->eval_diffusion_coefficients(*(Rght.gas));
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of one_d_interior_interp()

