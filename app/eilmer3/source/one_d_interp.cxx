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

// Some configuration settings...

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

bool apply_limiter = true;
bool set_apply_limiter_flag(bool bflag)
{
    apply_limiter = bflag;
    return apply_limiter;
}
bool get_apply_limiter_flag()
{
    return apply_limiter;
}

bool extrema_clipping = true;
bool set_extrema_clipping_flag(bool bflag)
{
    extrema_clipping = bflag;
    return extrema_clipping;
}
bool get_extrema_clipping_flag()
{
    return extrema_clipping;
}

// In the past, we've always done the high-order interpolation 
// of velocity in the global coordinate frame, however, the
// local frame is better.  Make that the default.
bool interpolate_in_local_frame = true;
bool set_interpolate_in_local_frame_flag(bool bflag)
{
    interpolate_in_local_frame = bflag;
    return interpolate_in_local_frame;
}
bool get_interpolate_in_local_frame_flag()
{
    return interpolate_in_local_frame;
}
   
//-----------------------------------------------------------------------------

/// \brief One-dimensional reconstruction of a scalar quantity.
///
/// See MBCNS workbook 2000/2 page 36 (26-Jan-2001) for formulation.
/// and MBCNS workbook 2005/Apr page 36 for new index labels

constexpr double epsilon = 1.0e-12; // Used within the van Albada limiter.

// The following module-level global variables must be set to appropriate values
// by one_d_interp_prepare() before their use in one_d_interp_scalar().
double aL0 = 0.0;
double aR0 = 0.0;
double lenL0_ = 0.0;
double lenR0_ = 0.0;
double two_over_lenL0_plus_lenL1 = 0.0;
double two_over_lenR0_plus_lenL0 = 0.0;
double two_over_lenR1_plus_lenR0 = 0.0;
double two_lenL0_plus_lenL1 = 0.0;
double two_lenR0_plus_lenR1 = 0.0;

inline int one_d_interp_both_prepare(double lenL1, double lenL0, double lenR0, double lenR1)
// Set up intermediate data that depends only on the cell geometry.
// It will remain constant when reconstructing the different scalar fields
// over the same set of cells.
{
    lenL0_ = lenL0;
    lenR0_ = lenR0;
    aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
    aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
    two_over_lenL0_plus_lenL1 = 2.0 / (lenL0 + lenL1);
    two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
    two_over_lenR1_plus_lenR0 = 2.0 / (lenR1 + lenR0);
    two_lenL0_plus_lenL1 = (2.0*lenL0 + lenL1);
    two_lenR0_plus_lenR1 = (2.0*lenR0 + lenR1);
    return SUCCESS;
} // end one_d_interp_both_prepare()

inline double clip_to_limits(double q, double A, double B)
// Returns q if q is between the values A and B, else
// it returns the closer limit of the range [A,B].
{
    double lower_limit = min(A, B);
    double upper_limit = max(A, B);
    return min(upper_limit, max(lower_limit, q));
} // end clip_to_limits()

inline int one_d_interp_both_scalar(double qL1, double qL0, double qR0, double qR1, double &qL, double &qR)
{
    double delLminus, del, delRplus, sL, sR;
    // Set up differences and limiter values.
    delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
    del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
    delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
    if ( apply_limiter ) {
	// val Albada limiter as per Ian Johnston's thesis.
	sL = (delLminus*del + fabs(delLminus*del)) / (delLminus*delLminus + del*del + epsilon);
	sR = (del*delRplus + fabs(del*delRplus)) / (del*del + delRplus*delRplus + epsilon);
    } else {
	// Use unlimited high-order reconstruction.
	sL = 1.0;
	sR = 1.0;
    }
    // The actual high-order reconstruction, possibly limited.
    qL = qL0 + sL * aL0 * ( del * two_lenL0_plus_lenL1 + delLminus * lenR0_ );
    qR = qR0 - sR * aR0 * ( delRplus * lenL0_ + del * two_lenR0_plus_lenR1 );
    if ( extrema_clipping ) {
	// An extra limiting filter to ensure that we do not compute new extreme values.
	// This was introduced to deal with very sharp transitions in species.
	qL = clip_to_limits(qL, qL0, qR0);
	qR = clip_to_limits(qR, qL0, qR0);
    }
    return SUCCESS;
} // end of one_d_interp_both_scalar()

inline int one_d_interp_left_prepare(double lenL1, double lenL0, double lenR0)
{
    lenL0_ = lenL0;
    lenR0_ = lenR0;
    aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
    two_over_lenL0_plus_lenL1 = 2.0 / (lenL0 + lenL1);
    two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
    two_lenL0_plus_lenL1 = (2.0*lenL0 + lenL1);
    return SUCCESS;
} // end one_d_interp_left_prepare()

inline int one_d_interp_left_scalar(double qL1, double qL0, double qR0, double &qL)
{
    double delLminus, del, sL;
    delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
    del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
    if ( apply_limiter ) {
	sL = (delLminus*del + fabs(delLminus*del)) / (delLminus*delLminus + del*del + epsilon);
    } else {
	sL = 1.0;
    }
    qL = qL0 + sL * aL0 * ( del * two_lenL0_plus_lenL1 + delLminus * lenR0_ );
    if ( extrema_clipping ) {
	qL = clip_to_limits(qL, qL0, qR0);
    }
    return SUCCESS;
} // end of one_d_interp_left_scalar()

inline int one_d_interp_right_prepare(double lenL0, double lenR0, double lenR1)
{
    lenL0_ = lenL0;
    lenR0_ = lenR0;
    aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
    two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
    two_over_lenR1_plus_lenR0 = 2.0 / (lenR1 + lenR0);
    two_lenR0_plus_lenR1 = (2.0*lenR0 + lenR1);
    return SUCCESS;
} // end one_d_interp_prepare()

inline int one_d_interp_right_scalar(double qL0, double qR0, double qR1, double &qR)
{
    double del, delRplus, sR;
    del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
    delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
    if ( apply_limiter ) {
	sR = (del*delRplus + fabs(del*delRplus)) / (del*del + delRplus*delRplus + epsilon);
    } else {
	sR = 1.0;
    }
    qR = qR0 - sR * aR0 * ( delRplus * lenL0_ + del * two_lenR0_plus_lenR1 );
    if ( extrema_clipping ) {
	qR = clip_to_limits(qR, qL0, qR0);
    }
    return SUCCESS;
} // end of one_d_interp_right_scalar()

//-----------------------------------------------------------------------------

/// \brief Reconstruct flow properties at an interface from a full set of 4 FV_Cell properties.
///
/// This is essentially a one-dimensional interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.
int one_d_interp_both(const FV_Interface &IFace,
		      const FV_Cell &cL1, const FV_Cell &cL0,
		      const FV_Cell &cR0, const FV_Cell &cR1,
		      double cL1Length, double cL0Length,
		      double cR0Length, double cR1Length,
		      FlowState &Lft, FlowState &Rght )
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();

    // Low-order reconstruction just copies data from adjacent FV_Cell.
    // Even for high-order reconstruction, we depend upon this copy for
    // the viscous-transport and diffusion coefficients.
    Lft.copy_values_from(*(cL0.fs));
    Rght.copy_values_from(*(cR0.fs));
    if ( G.Xorder > 1 ) {
	// High-order reconstruction for some properties.
	if ( interpolate_in_local_frame ) {
	    // Paul Petrie-Repar and Jason Qin have noted that the velocity needs
	    // to be reconstructed in the interface-local frame of reference so that
	    // the normal velocities are not messed up for mirror-image at walls.
	    // PJ 21-feb-2012
	    cL1.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	    cL0.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	    cR0.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	    cR1.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	}
	one_d_interp_both_prepare(cL1Length, cL0Length, cR0Length, cR1Length);
	one_d_interp_both_scalar(cL1.fs->vel.x, cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x, Lft.vel.x, Rght.vel.x);
	one_d_interp_both_scalar(cL1.fs->vel.y, cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y, Lft.vel.y, Rght.vel.y);
	one_d_interp_both_scalar(cL1.fs->vel.z, cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z, Lft.vel.z, Rght.vel.z);
	if ( G.MHD ) {
	    one_d_interp_both_scalar(cL1.fs->B.x, cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x, Lft.B.x, Rght.B.x);
	    one_d_interp_both_scalar(cL1.fs->B.y, cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y, Lft.B.y, Rght.B.y);
	    one_d_interp_both_scalar(cL1.fs->B.z, cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z, Lft.B.z, Rght.B.z);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    one_d_interp_both_scalar(cL1.fs->tke, cL0.fs->tke, cR0.fs->tke, cR1.fs->tke, Lft.tke, Rght.tke);
	    one_d_interp_both_scalar(cL1.fs->omega, cL0.fs->omega, cR0.fs->omega, cR1.fs->omega, Lft.omega, Rght.omega);
	}
	Gas_data &gL1 = *(cL1.fs->gas); Gas_data &gL0 = *(cL0.fs->gas);
	Gas_data &gR0 = *(cR0.fs->gas); Gas_data &gR1 = *(cR1.fs->gas);
	if ( nsp > 1 ) {
	    // Multiple species.
	    for ( size_t isp = 0; isp < nsp; ++isp ) {
		one_d_interp_both_scalar(gL1.massf[isp], gL0.massf[isp], gR0.massf[isp], gR1.massf[isp],
					 Lft.gas->massf[isp], Rght.gas->massf[isp]);
	    }
	    if ( scale_mass_fractions( Lft.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp < nsp; ++isp )
		    Lft.gas->massf[isp] = gL0.massf[isp];
	    }
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp < nsp; ++isp )
		    Rght.gas->massf[isp] = gR0.massf[isp];
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
	    one_d_interp_both_scalar(gL1.p, gL0.p, gR0.p, gR1.p, Lft.gas->p, Rght.gas->p);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_both_scalar(gL1.T[i], gL0.T[i], gR0.T[i], gR1.T[i], Lft.gas->T[i], Rght.gas->T[i]);
	    }
	    if ( gmodel->eval_thermo_state_pT(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    if ( gmodel->eval_thermo_state_pT(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOE:
	    one_d_interp_both_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas->rho, Rght.gas->rho);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_both_scalar(gL1.e[i], gL0.e[i], gR0.e[i], gR1.e[i], Lft.gas->e[i], Rght.gas->e[i]);
	    }
	    if ( gmodel->eval_thermo_state_rhoe(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    if ( gmodel->eval_thermo_state_rhoe(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOP:
	    one_d_interp_both_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas->rho, Rght.gas->rho);
	    one_d_interp_both_scalar(gL1.p, gL0.p, gR0.p, gR1.p, Lft.gas->p, Rght.gas->p);
	    if ( gmodel->eval_thermo_state_rhop(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    if ( gmodel->eval_thermo_state_rhop(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOT: 
	    one_d_interp_both_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas->rho, Rght.gas->rho);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_both_scalar(gL1.T[i], gL0.T[i], gR0.T[i], gR1.T[i], Lft.gas->T[i], Rght.gas->T[i]);
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
	if ( interpolate_in_local_frame ) {
	    // Undo the transformation made earlier. PJ 21-feb-2012
	    Lft.vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    Rght.vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cL1.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cL0.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cR0.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cR1.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of one_d_interp_both()

/// \brief Reconstruct flow properties at an interface from a cells L1,L0,R0.
///
/// This is essentially a one-dimensional interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.
int one_d_interp_left(const FV_Interface &IFace,
		      const FV_Cell &cL1, const FV_Cell &cL0, const FV_Cell &cR0,
		      double cL1Length, double cL0Length, double cR0Length,
		      FlowState &Lft, FlowState &Rght)
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();

    // Low-order reconstruction just copies data from adjacent FV_Cell.
    // Even for high-order reconstruction, we depend upon this copy for
    // the viscous-transport and diffusion coefficients.
    Lft.copy_values_from(*(cL0.fs));
    Rght.copy_values_from(*(cR0.fs));
    if ( G.Xorder > 1 ) {
	// High-order reconstruction for some properties.
	if ( interpolate_in_local_frame ) {
	    // In the interface-local frame.
	    cL1.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	    cL0.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	    cR0.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	}
	one_d_interp_left_prepare(cL1Length, cL0Length, cR0Length);
	one_d_interp_left_scalar(cL1.fs->vel.x, cL0.fs->vel.x, cR0.fs->vel.x, Lft.vel.x);
	one_d_interp_left_scalar(cL1.fs->vel.y, cL0.fs->vel.y, cR0.fs->vel.y, Lft.vel.y);
	one_d_interp_left_scalar(cL1.fs->vel.z, cL0.fs->vel.z, cR0.fs->vel.z, Lft.vel.z);
	if ( G.MHD ) {
	    one_d_interp_left_scalar(cL1.fs->B.x, cL0.fs->B.x, cR0.fs->B.x, Lft.B.x);
	    one_d_interp_left_scalar(cL1.fs->B.y, cL0.fs->B.y, cR0.fs->B.y, Lft.B.y);
	    one_d_interp_left_scalar(cL1.fs->B.z, cL0.fs->B.z, cR0.fs->B.z, Lft.B.z);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    one_d_interp_left_scalar(cL1.fs->tke, cL0.fs->tke, cR0.fs->tke, Lft.tke);
	    one_d_interp_left_scalar(cL1.fs->omega, cL0.fs->omega, cR0.fs->omega, Lft.omega);
	}
	Gas_data &gL1 = *(cL1.fs->gas); Gas_data &gL0 = *(cL0.fs->gas); Gas_data &gR0 = *(cR0.fs->gas);
	if ( nsp > 1 ) {
	    // Multiple species.
	    for ( size_t isp = 0; isp < nsp; ++isp ) {
		one_d_interp_left_scalar(gL1.massf[isp], gL0.massf[isp], gR0.massf[isp], Lft.gas->massf[isp]);
	    }
	    if ( scale_mass_fractions( Lft.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp < nsp; ++isp )
		    Lft.gas->massf[isp] = gL0.massf[isp];
	    }
	} else {
	    // Only one possible mass-fraction value for a single species.
	    Lft.gas->massf[0] = 1.0;
	}
	// Interpolate on two of the thermodynamic quantities, 
	// and fill in the rest based on an EOS call. 
	// If an EOS call fails, fall back to just copying cell-centre data.
	// This does presume that the cell-centre data is valid. 
	switch ( thermo_interpolator ) {
	case INTERP_PT: 
	    one_d_interp_left_scalar(gL1.p, gL0.p, gR0.p, Lft.gas->p);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_left_scalar(gL1.T[i], gL0.T[i], gR0.T[i], Lft.gas->T[i]);
	    }
	    if ( gmodel->eval_thermo_state_pT(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    break;
	case INTERP_RHOE:
	    one_d_interp_left_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas->rho);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_left_scalar(gL1.e[i], gL0.e[i], gR0.e[i], Lft.gas->e[i]);
	    }
	    if ( gmodel->eval_thermo_state_rhoe(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    break;
	case INTERP_RHOP:
	    one_d_interp_left_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas->rho);
	    one_d_interp_left_scalar(gL1.p, gL0.p, gR0.p, Lft.gas->p);
	    if ( gmodel->eval_thermo_state_rhop(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    break;
	case INTERP_RHOT: 
	    one_d_interp_left_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas->rho);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_left_scalar(gL1.T[i], gL0.T[i], gR0.T[i], Lft.gas->T[i]);
	    }
	    if ( gmodel->eval_thermo_state_rhoT(*(Lft.gas)) != SUCCESS ) {
		Lft.copy_values_from(*(cL0.fs));
	    }
	    break;
	default: 
	    throw std::runtime_error("Invalid thermo interpolator.");
	}
	if ( interpolate_in_local_frame ) {
	    // Undo the transformation made earlier.
	    Lft.vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cL1.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cL0.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cR0.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of one_d_interp_left()

/// \brief Reconstruct flow properties at an interface from cells L0,R0,R1.
///
/// This is essentially a one-dimensional interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.
int one_d_interp_right(const FV_Interface &IFace,
		       const FV_Cell &cL0, const FV_Cell &cR0, const FV_Cell &cR1,
		       double cL0Length, double cR0Length, double cR1Length,
		       FlowState &Lft, FlowState &Rght)
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();

    // Low-order reconstruction just copies data from adjacent FV_Cell.
    // Even for high-order reconstruction, we depend upon this copy for
    // the viscous-transport and diffusion coefficients.
    Lft.copy_values_from(*(cL0.fs));
    Rght.copy_values_from(*(cR0.fs));
    if ( G.Xorder > 1 ) {
	// High-order reconstruction for some properties.
	if ( interpolate_in_local_frame ) {
	    // In the interface-local frame.
	    cL0.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	    cR0.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	    cR1.fs->vel.transform_to_local(IFace.n, IFace.t1, IFace.t2);
	}
	one_d_interp_right_prepare(cL0Length, cR0Length, cR1Length);
	one_d_interp_right_scalar(cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x, Rght.vel.x);
	one_d_interp_right_scalar(cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y, Rght.vel.y);
	one_d_interp_right_scalar(cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z, Rght.vel.z);
	if ( G.MHD ) {
	    one_d_interp_right_scalar(cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x, Rght.B.x);
	    one_d_interp_right_scalar(cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y, Rght.B.y);
	    one_d_interp_right_scalar(cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z, Rght.B.z);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    one_d_interp_right_scalar(cL0.fs->tke, cR0.fs->tke, cR1.fs->tke, Rght.tke);
	    one_d_interp_right_scalar(cL0.fs->omega, cR0.fs->omega, cR1.fs->omega, Rght.omega);
	}
	Gas_data &gL0 = *(cL0.fs->gas); Gas_data &gR0 = *(cR0.fs->gas); Gas_data &gR1 = *(cR1.fs->gas);
	if ( nsp > 1 ) {
	    // Multiple species.
	    for ( size_t isp = 0; isp < nsp; ++isp ) {
		one_d_interp_right_scalar(gL0.massf[isp], gR0.massf[isp], gR1.massf[isp], Rght.gas->massf[isp]);
	    }
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp < nsp; ++isp )
		    Rght.gas->massf[isp] = gR0.massf[isp];
	    }
	} else {
	    // Only one possible mass-fraction value for a single species.
	    Rght.gas->massf[0] = 1.0;
	}
	// Interpolate on two of the thermodynamic quantities, 
	// and fill in the rest based on an EOS call. 
	// If an EOS call fails, fall back to just copying cell-centre data.
	// This does presume that the cell-centre data is valid. 
	switch ( thermo_interpolator ) {
	case INTERP_PT: 
	    one_d_interp_right_scalar(gL0.p, gR0.p, gR1.p, Rght.gas->p);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_right_scalar(gL0.T[i], gR0.T[i], gR1.T[i], Rght.gas->T[i]);
	    }
	    if ( gmodel->eval_thermo_state_pT(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOE:
	    one_d_interp_right_scalar(gL0.rho, gR0.rho, gR1.rho, Rght.gas->rho);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_right_scalar(gL0.e[i], gR0.e[i], gR1.e[i], Rght.gas->e[i]);
	    }
	    if ( gmodel->eval_thermo_state_rhoe(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOP:
	    one_d_interp_right_scalar(gL0.rho, gR0.rho, gR1.rho, Rght.gas->rho);
	    one_d_interp_right_scalar(gL0.p, gR0.p, gR1.p, Rght.gas->p);
	    if ( gmodel->eval_thermo_state_rhop(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	case INTERP_RHOT: 
	    one_d_interp_right_scalar(gL0.rho, gR0.rho, gR1.rho, Rght.gas->rho);
	    for ( size_t i = 0; i < nmodes; ++i ) {
		one_d_interp_right_scalar(gL0.T[i], gR0.T[i], gR1.T[i], Rght.gas->T[i]);
	    }
	    if ( gmodel->eval_thermo_state_rhoT(*(Rght.gas)) != SUCCESS ) {
		Rght.copy_values_from(*(cR0.fs));
	    }
	    break;
	default: 
	    throw std::runtime_error("Invalid thermo interpolator.");
	}
	if ( interpolate_in_local_frame ) {
	    // Undo the transformation made earlier.
	    Rght.vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cL0.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cR0.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	    cR1.fs->vel.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of one_d_interp_right()

