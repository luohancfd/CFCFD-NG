/** \file LOS_pieces.cxx
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 09-Dec-09
 *
 **/

#include <stdlib.h> 
#include <math.h>
#include <sstream>
#include <iostream>

#include "../../nm/source/exponential_integrals.hh"
#include "../../util/source/useful.h"

#include "LOS_pieces.hh"

using namespace std;

RadiatingPoint::RadiatingPoint()
{
    // 1. Initialise the CoeffSpectra
    X_ = new CoeffSpectra();
    
    // gas-data structure etc will be defined later
}

RadiatingPoint::RadiatingPoint( Gas_data * Q, double * Q_rE_rad, double s, double ds )
{
    // just call the redefine function
    this->redefine(Q,Q_rE_rad,s,ds);
}

RadiatingPoint::~RadiatingPoint()
{
    delete X_;
}

void RadiatingPoint::redefine( Gas_data * Q, double * Q_rE_rad, double s, double ds )
{
    // 1. Set pointers and s val
    Q_ = Q;
    Q_rE_rad_ = Q_rE_rad;
    s_ = s;
    ds_ = ds;
    
    return;
}

/* ---------- Line-of-sight data storage and manipulation ----------- */

LOS_data::LOS_data( RadiationSpectralModel * rsm, int nrps, double T_i, double T_f )
 : rsm_( rsm ), nrps_( nrps ), T_i_( T_i ), T_f_( T_f )
{
    rpoints_.resize( nrps_ );
    for ( int irp=0; irp<nrps_; ++irp )
    	rpoints_[irp] = new RadiatingPoint();
    // Assuming uniform spectral distribution
    nnus_ = rsm_->get_spectral_points();
}
 
LOS_data::~LOS_data()
{
    for ( size_t irp=0; irp<rpoints_.size(); ++irp )
    	delete rpoints_[irp];
}

void LOS_data::clear()
{
    for ( size_t irp=0; irp<rpoints_.size(); ++irp )
    	delete rpoints_[irp];
}

void LOS_data::set_rad_point( int irp, Gas_data * Q, double * Q_rE_rad, double s, double ds )
{
    rpoints_[irp]->redefine( Q, Q_rE_rad, s, ds );
    rsm_->radiative_spectra_for_gas_state( *(rpoints_[irp]->Q_), *(rpoints_[irp]->X_) );
}

void LOS_data::clone_rad_point( int iprp, int irp, double * Q_rE_rad, double s, double ds )
{
    rpoints_[irp]->redefine( rpoints_[iprp]->Q_, Q_rE_rad, s, ds );
    delete rpoints_[irp]->X_;
    rpoints_[irp]->X_ = rpoints_[iprp]->X_->clone();
}

double LOS_data::integrate_LOS( SpectralIntensity &S )
{
    /* Spectral LOS integration using the exact transport eq. sol. for constant property slabs */
    
    // Check vector sizes in S for consistency
    if ( int(S.nu.size()) != nnus_ || int(S.I_nu.size()) != nnus_ ) {
	cout << "Mismatch in size of spectral vector in structure S. " << endl
	     << "nnus = " << nnus_ << ", S.nu.size() = " << S.nu.size() << endl
	     << "Bailing out!" << endl;
	exit( BAD_INPUT_ERROR );
    }
    
    // Compute and append Planck backlight intensity
    for ( int inu=0; inu < nnus_; inu++ )
	S.I_nu[inu] += planck_intensity( S.nu[inu], T_i_ );
    
    // Spatial integration over discretised points
    for ( int irp=0; irp<nrps_; irp++ ) {
	double ds = rpoints_[irp]->ds_;
	// Spatial integration for each spectral interval
	for ( int inu=0; inu < nnus_; inu++ ) {
	    // exact solution for constant property slab
	    double decay_factor = exp( - ds * rpoints_[irp]->X_->kappa_nu[inu] );
	    // newly radiated intensity
	    double tmpA = rpoints_[irp]->X_->j_nu[inu] / rpoints_[irp]->X_->kappa_nu[inu] * ( 1.0 - decay_factor );
	    if ( isnan(tmpA) ) {
		// tmpA will go NaN when kappa_nu is zero.  Intensity is zero.
		tmpA = 0.0;
	    }
	    // incoming intensity
	    double tmpB = S.I_nu[inu] * decay_factor;
	    // write-over old I_nu with new I_nu
	    S.I_nu[inu] = tmpA + tmpB;
	}
    }
    // Subtract out Planck intensity from outer BC and cumpute integrated intensity
    double I_total = 0.0;
    for ( int inu=0; inu < nnus_; inu++ ) {
	S.I_nu[inu] -= planck_intensity( S.nu[inu], T_f_ );
	if ( inu>0 )
	    I_total += 0.5 * ( S.I_nu[inu] + S.I_nu[inu-1] ) * ( S.nu[inu] - S.nu[inu-1] );
    }

    // End void function
    return I_total;
}

void LOS_data::write_all_points_to_file( void )
{
    for ( size_t ip=0; ip<rpoints_.size(); ++ip ) {
    	ostringstream fname;
    	fname << "rpoint-" << ip << ".dat";
    	rpoints_[ip]->X_->write_to_file( fname.str() );
    }
    
    return;
}

void LOS_data::write_point_to_file( int ip, string fname )
{
    if ( ip >= nrps_ ) {
    	cout << "LOS_data::write_points_to_file()" << endl
    	     << "ip = " << ip << ", nrps_ = " << nrps_ << endl
    	     << "Bailing out!" << endl;
    	exit( BAD_INPUT_ERROR );
    }
    rpoints_[ip]->X_->write_to_file( fname );
    
    return;
}

/* ---------- Tangent-slab data storage and manipulation ----------- */

TS_data::TS_data( RadiationSpectralModel * rsm, int nrps, double T_i, double T_f )
 : LOS_data( rsm, nrps, T_i, T_f )
{
    F_ = new SpectralFlux( rsm_ );
    F_->q_nu.resize( nnus_ );
    F_->nu.resize( nnus_ );
    /* Uniformally distributed spectral points with constant frequency spacing */
    double nu = lambda2nu( rsm_->get_lambda_max() );
    double dnu = ( lambda2nu( rsm_->get_lambda_min() ) - lambda2nu( rsm_->get_lambda_max() ) ) 
		/ double ( ( rsm_->get_spectral_points() - 1 ) );
    F_->nu.resize( rsm_->get_spectral_points() );
    for( int inu=0; inu<rsm_->get_spectral_points(); ++inu ){
	F_->nu[inu] = nu;
	nu+=dnu;
    }
}
 
TS_data::~TS_data()
{
    delete F_;
}

void TS_data::clear()
{
    for ( size_t inu=0; inu<F_->nu.size(); ++inu )
    	F_->q_nu[inu] = 0.0;
    
    for ( size_t irp=0; irp<rpoints_.size(); ++irp )
    	delete rpoints_[irp];
}

double TS_data::quick_solve_for_divq()
{
    // Simple method from Sebastian Karl's VKI thesis where q+ and q- between 
    // each node location are solved for
    
    // Make vectors to store integrated fluxes
    vector<double> q_plus, q_minus;
    
    /* ------- Wall-to-shock (minus) direction ------- */
    
    // Outer-boundary blackbody intensity (note optical thickness = 0)
    double q_int = 0.0;
    for ( int inu=0; inu < nnus_; inu++ ) {
        F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_f_ ) * E_3( 0.0 );
        // Frequency integration
        if ( inu>0 )
            q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
    }
    
    // Save integrated flux
    q_minus.push_back( q_int );
    
    // integrate over emitting slabs using trapezoidal method
    for ( int irp=(nrps_-1); irp >= 0; --irp ) {
        // cout << "Spatial location: | " << irp << " | \r" << flush;
        q_int = 0.0;
        for ( int inu=0; inu < nnus_; inu++ ) {
            // emission and absorption coefficients for this slab
            double kappa_nu = rpoints_[irp]->X_->kappa_nu[inu];
            double j_nu = rpoints_[irp]->X_->j_nu[inu];
            // calculate optical thickness between intervals
            double tau_nu = fabs(rpoints_[irp]->ds_) * kappa_nu;
            // Attenuate incident flux )
            F_->q_nu[inu] *= 2.0 * E_3( tau_nu );
            // Add contribution from this slab
	    double tmpA = 0.0;
	    if ( j_nu != 0.0 )	tmpA = j_nu / kappa_nu;
	     F_->q_nu[inu] += tmpA * M_PI * ( 1.0 -  2.0 * E_3( tau_nu ) );
            // Integrate
            if ( inu>0 )
                q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
        }
        // Save integrated flux
        q_minus.push_back( q_int );
    }
    
    /* ------- Shock-to-wall (plus) direction ------- */
    
    // Outer-boundary blackbody intensity (note optical thickness = 0)
    q_int = 0.0;
    for ( int inu=0; inu < nnus_; inu++ ) {
        F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_i_ ) * E_3( 0.0 );
        // Frequency integration
        if ( inu>0 )
            q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
    }
    
    // Save integrated flux
    q_plus.push_back( q_int );
    
    // integrate over emitting slabs using trapezoidal method
    for ( int irp=0; irp < nrps_; irp++ ) {
        // cout << "Spatial location: | " << irp << " | \r" << flush;
        q_int = 0.0;
        for ( int inu=0; inu < nnus_; inu++ ) {
            // emission and absorption coefficients for this slab
            double kappa_nu = rpoints_[irp]->X_->kappa_nu[inu];
            double j_nu = rpoints_[irp]->X_->j_nu[inu];
            // calculate optical thickness between intervals
            double tau_nu = fabs(rpoints_[irp]->ds_) * kappa_nu;
            // Attenuate incident flux )
            F_->q_nu[inu] *= 2.0 * E_3( tau_nu );
            // Add contribution from this slab
	    double tmpA = 0.0;
	    if ( j_nu != 0.0 )	tmpA = j_nu / kappa_nu;
	     F_->q_nu[inu] += tmpA * M_PI * ( 1.0 -  2.0 * E_3( tau_nu ) );
            // Integrate
            if ( inu>0 )
                q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
        }
        // Save integrated flux
        q_plus.push_back( q_int );
    }

    // calculate radiative divergence

    for ( int irp=0; irp<nrps_; ++irp ) {
	// use forward differencing to evaluate divergence
	double q_net_i = q_plus[irp] - q_minus[nrps_-irp];
	double q_net_ip1 = q_plus[irp+1] - q_minus[nrps_-irp-1];
	*(rpoints_[irp]->Q_rE_rad_) = ( q_net_i - q_net_ip1 ) / rpoints_[irp]->ds_;
    }
    
    // return wall directed flux
    // NOTE: omitting wall emission - which is correct
    double q_wall = q_plus.back();
    
    return q_wall;
}

double TS_data::exact_solve_for_divq()
{
    /* Exact tangent-slab equations as presented by Modest and Feldick et al */
    
    // Make a vector to store integrated fluxes
    vector<double> q_plus, q_minus;
    
    /* ------- Wall-to-shock (minus) direction ------- */
    
    // Outer-boundary blackbody intensity (note optical thickness = 0)
    double q_int = 0.0;
    for ( int inu=0; inu < nnus_; inu++ ) {
        F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_f_ ) * E_3( 0.0 );
        // Frequency integration
        if ( inu>0 )
            q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
    }
    
    // save flux at wall from wall
    q_minus.push_back( q_int );
    
    for ( int krp=(nrps_-1); krp >= 0; --krp ) {
	q_int = 0.0;
	for ( int inu=0; inu < nnus_; inu++ ) {
	    double ds_i = rpoints_[nrps_-1]->ds_;
	    // Firstly do contribution from wall
	    double dtau_k = 0.0;
	    double ds_j = ds_i;
	    for ( int jrp=krp; jrp >= 0; --jrp ) {
		dtau_k += fabs(ds_j) * rpoints_[jrp]->X_->kappa_nu[inu];
		// update ds_j
		ds_j = rpoints_[jrp]->ds_;
	    }
	    F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_f_ ) * E_3( dtau_k );
	    // now contributions from cells, starting from closest
	    for ( int irp=(nrps_-1); irp >= krp; --irp ) {
		// source function for original cell
		double S = 0.0;
		if ( rpoints_[irp]->X_->j_nu[inu] != 0.0 ) S = rpoints_[irp]->X_->j_nu[inu] / rpoints_[irp]->X_->kappa_nu[inu];
		double dtau_jm1 = 0.0;
		ds_j = ds_i;
		for ( int jrp=irp; jrp >= krp; --jrp ) {
		    dtau_jm1 += fabs(ds_j) * rpoints_[jrp]->X_->kappa_nu[inu];
		    // update ds_j
		    ds_j = rpoints_[jrp]->ds_;
		}
		double dtau_j = dtau_jm1 - rpoints_[irp]->X_->kappa_nu[inu] * fabs(ds_i);
		// add contribution from this cell
		F_->q_nu[inu] += 2.0 * M_PI * S * ( E_3( dtau_j ) -  E_3( dtau_jm1 ) );
		// update ds_i
		if ( irp != 0 )
		    ds_i = 2.0 * ( rpoints_[irp-1]->s_ - rpoints_[irp]->s_ - 0.5 * ds_i );
	    }
            if ( inu>0 )
                q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
	}
	q_minus.push_back( q_int );
    }
    
    /* ------- Shock-to-wall (plus) direction ------- */
    
    // Outer-boundary blackbody intensity (note optical thickness = 0)
    q_int = 0.0;
    for ( int inu=0; inu < nnus_; inu++ ) {
        F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_i_ ) * E_3( 0.0 );
        // Frequency integration
        if ( inu>0 )
            q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
    }
    
    // Save integrated flux
    q_plus.push_back( q_int );
    
    for ( int krp=0; krp < nrps_; ++krp ) {
	q_int = 0.0;
	for ( int inu=0; inu < nnus_; inu++ ) {
	    double ds_i = rpoints_[0]->ds_;
	    // Firstly do contribution from wall
	    double dtau_k = 0.0;
	    double ds_j = ds_i;
	    for ( int jrp=0; jrp <= krp; ++jrp ) {
		dtau_k += fabs(ds_j) * rpoints_[jrp]->X_->kappa_nu[inu];
		// update ds_j
		ds_j = 2.0 * rpoints_[jrp]->ds_;
	    }
	    F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_i_ ) * E_3( dtau_k );
	    // now contributions from cells, starting from closest
	    for ( int irp=0; irp <= krp; ++irp ) {
		// source function for original cell
		double S = 0.0;
		if ( rpoints_[irp]->X_->j_nu[inu] != 0.0 ) S = rpoints_[irp]->X_->j_nu[inu] / rpoints_[irp]->X_->kappa_nu[inu];
		double dtau_jm1 = 0.0;
		double ds_j = ds_i;
		for ( int jrp=irp; jrp <= krp; jrp++ ) {
		    dtau_jm1 += ds_j * rpoints_[jrp]->X_->kappa_nu[inu];
		    // update ds_j
		    ds_j = rpoints_[jrp]->ds_;
		}
		double dtau_j = dtau_jm1 - rpoints_[irp]->X_->kappa_nu[inu] * ds_i;
		// add contribution from this cell
		F_->q_nu[inu] += 2.0 * M_PI * S * ( E_3( dtau_j ) -  E_3( dtau_jm1 ) );
		// update ds_i
		ds_i = rpoints_[irp]->ds_;
	    }
            if ( inu>0 )
                q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
	}
	q_plus.push_back( q_int );
    }
    
    // calculate radiative divergence
    for ( int irp=0; irp<nrps_; ++irp ) {
    	double ds = rpoints_[irp]->ds_;
	// use forward differencing to evaluate divergence
	double q_net_i = q_plus[irp] - q_minus[nrps_-irp];
	double q_net_ip1 = q_plus[irp+1] - q_minus[nrps_-irp-1];
	*(rpoints_[irp]->Q_rE_rad_) = ( q_net_i - q_net_ip1 ) / ds;
	// cout << "q_plus = " << q_plus[irp] << ", q_minus = " << q_minus[irp] << endl;
    }
    
    // return wall directed flux
    // NOTE: omitting wall emission - which is correct
    double q_wall = q_plus.back();
    
    return q_wall;
}

