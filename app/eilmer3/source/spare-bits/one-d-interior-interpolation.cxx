// This code pulled out of one_d_interp.cxx 24-Aug-2013
// It was adapted from the other one_d_interp functions by Andrew Pastrello
// for his shock-fitting work.


/// \brief Asymmetric one-dimensional reconstruction of a scalar quantity.
///
/// See Ian Johnston's thesis.
///
inline int one_d_interior_interp_scalar( double qL0, double qR0, double qR1, double qR2, 
					 double lenL0, double lenR0, double lenR1, double lenR2,
					 double &qL, double &qR )
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
    if ( apply_limiter ) {
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
    if ( extrema_clipping ) {
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

    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Lft.copy_values_from(*(cL0.fs));
    Rght.copy_values_from(*(cR0.fs));
    if ( G.Xorder > 1 ) {
	// High-order reconstruction for some properties.
	one_d_interior_interp_scalar(cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x, cR2.fs->vel.x,
				     cL0Length, cR0Length, cR1Length, cR2Length,
				     Lft.vel.x, Rght.vel.x);
	one_d_interior_interp_scalar(cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y, cR2.fs->vel.y,
				     cL0Length, cR0Length, cR1Length, cR2Length,
				     Lft.vel.y, Rght.vel.y);
	one_d_interior_interp_scalar(cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z, cR2.fs->vel.z,
				     cL0Length, cR0Length, cR1Length, cR2Length,
				     Lft.vel.z, Rght.vel.z);
	if ( get_mhd_flag() == 1 ) {
	    one_d_interior_interp_scalar(cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x, cR2.fs->B.x,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.B.x, Rght.B.x);
	    one_d_interior_interp_scalar(cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y, cR2.fs->B.y,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.B.y, Rght.B.y);
	    one_d_interior_interp_scalar(cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z, cR2.fs->B.z,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.B.z, Rght.B.z);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    one_d_interior_interp_scalar(cL0.fs->tke, cR0.fs->tke, cR1.fs->tke, cR2.fs->tke,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.tke, Rght.tke);
	    one_d_interior_interp_scalar(cL0.fs->omega, cR0.fs->omega, cR1.fs->omega, cR2.fs->omega,
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.omega, Rght.omega);
	}
        for ( int isp = 0; isp < nsp; ++isp ) {
	    one_d_interior_interp_scalar(cL0.fs->gas->massf[isp], cR0.fs->gas->massf[isp],
					 cR1.fs->gas->massf[isp], cR2.fs->gas->massf[isp],
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.gas->massf[isp], Rght.gas->massf[isp]);
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
				     Lft.gas->rho, Rght.gas->rho);

	for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
	    one_d_interior_interp_scalar(gL0.e[i], gR0.e[i], gR1.e[i], gR2.e[i],
					 cL0Length, cR0Length, cR1Length, cR2Length,
					 Lft.gas->e[i], Rght.gas->e[i]);
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of one_d_interior_interp()

