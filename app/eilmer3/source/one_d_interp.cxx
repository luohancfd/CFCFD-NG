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
#include "one_d_interp_scalar.hh"
#include "cell.hh"
#include "one_d_interp.hh"
#include "kernel.hh"
// #include "flux_calc.hh"
// #include "diffusion.hh"
// #include "bgk.hh"

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
    Thermo_interpolator* ti = get_thermo_interpolator_ptr();
    int nsp = gmodel->get_number_of_species();
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
        for ( int isp = 0; isp < nsp; ++isp ) {
	    one_d_interp_scalar(cL1.fs->gas->massf[isp], cL0.fs->gas->massf[isp],
				cR0.fs->gas->massf[isp], cR1.fs->gas->massf[isp],
				cL1Length, cL0Length, cR0Length, cR1Length,
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

	int status = ti->one_d_interp(*(cL1.fs->gas), *(cL0.fs->gas), *(cR0.fs->gas), *(cR1.fs->gas),
				      cL1Length, cL0Length, cR0Length, cR1Length,
				      *(Lft.gas), *(Rght.gas));
	if ( status != SUCCESS ) {
	    if ( status == 1 ) {
		// Lft state failed.
		Lft.copy_values_from(*(cL0.fs));
	    }
	    else if ( status == 2 ) {
		// Rght state failed.
		Rght.copy_values_from(*(cR0.fs));
	    }
	    else if ( status == 3 ) {
		// Both failed.
		Lft.copy_values_from(*(cL0.fs));
		Rght.copy_values_from(*(cR0.fs));
	    }
	    else {
		printf("one_d_interp(): Problem in flow state reconstruction.");
		printf("Failure status: %d is unknown.\n", status);
		printf("Bailing out!\n");
		exit(RECONSTRUCTION_ERROR);
	    }
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
} // end of one_d_interp()

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
	int status = gmodel->eval_thermo_state_rhoe(*(Lft.gas));
	status = gmodel->eval_thermo_state_rhoe(*(Rght.gas));
	if ( status != SUCCESS ) {
	    if ( status == 1 ) {
		// Lft state failed.
		Lft.copy_values_from(*(cL0.fs));
	    }
	    else if ( status == 2 ) {
		// Rght state failed.
		Rght.copy_values_from(*(cR0.fs));
	    }
	    else if ( status == 3 ) {
		// Both failed.
		Lft.copy_values_from(*(cL0.fs));
		Rght.copy_values_from(*(cR0.fs));
	    }
	    else {
		printf("one_d_interp(): Problem in flow state reconstruction.");
		printf("Failure status: %d is unknown.\n", status);
		printf("Bailing out!\n");
		exit(RECONSTRUCTION_ERROR);
	    }
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

