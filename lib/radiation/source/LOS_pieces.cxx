/** \file LOS_pieces.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 09-Dec-09
 *
 **/

#include <stdlib.h> 
#include <math.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../nm/source/exponential_integrals.hh"
#include "../../util/source/useful.h"
#include "../../util/source/randomc.h"

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

LOS_data::LOS_data( RadiationSpectralModel * rsm, int nrps, double T_i, double T_f, bool store_spectra_on_disk )
 : rsm_( rsm ), nrps_( nrps ), T_i_( T_i ), T_f_( T_f ), store_spectra_on_disk_( store_spectra_on_disk )
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
    if ( store_spectra_on_disk_ ) {
        ostringstream oss;
        oss << "rpoint_" << irp << ".txt";
        rpoints_[irp]->X_->write_to_file(oss.str(),FREQUENCY);
        rpoints_[irp]->X_->clear_data();
    }
}

void LOS_data::clone_rad_point( int iprp, int irp, double * Q_rE_rad, double s, double ds )
{
    rpoints_[irp]->redefine( rpoints_[iprp]->Q_, Q_rE_rad, s, ds );
    delete rpoints_[irp]->X_;
    rpoints_[irp]->X_ = rpoints_[iprp]->X_->clone();
}

void LOS_data::create_spectral_bins( int binning_type, int N_bins, vector<SpectralBin*> & B )
{
    // check that we have some rpoints_
    if ( rpoints_.size()==0 ) {
	cout << "LOS_data::create_spectral_bins()" << endl
	     << "No spectral points have been created." << endl
	     << "Exiting program." << endl;
	exit(FAILURE);
    }

    if ( binning_type==FREQUENCY_BINNING ) {
	create_spectral_bin_vector( rpoints_[0]->X_->nu, binning_type, N_bins, B );
    }
    else if ( binning_type==OPACITY_BINNING ) {
	// We need to solve for the spatially independent mean opacity/absorption (see Eq 2.2 of Wray, Ripoll and Prabhu)
	// NOTE: We are assuming planar geometry (ie that the cells have the same width in the z direction)
	//       For axi geometries, this is only accurate along stagnation streamlines
	vector<double> kappa_mean;
	for ( int inu=0; inu < nnus_; inu++ ) {
	    double j_nu_dV = 0.0, S_nu_dV = 0.0;
	    for ( int irp=0; irp<nrps_; irp++ ) {
		double dj_nu_dV = rpoints_[irp]->X_->j_nu[inu] * rpoints_[irp]->ds_;
		j_nu_dV += dj_nu_dV;
		if ( rpoints_[irp]->X_->kappa_nu[inu] > 0.0 )
		    S_nu_dV += dj_nu_dV / rpoints_[irp]->X_->kappa_nu[inu];
	    }
	    // If the source function is zero then kappa should also be zero
	    if ( S_nu_dV==0.0 )
		kappa_mean.push_back( 0.0 );
	    else
		kappa_mean.push_back( j_nu_dV / S_nu_dV );
	}
	create_spectral_bin_vector( kappa_mean, binning_type, N_bins, B );
    }
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
    if ( nnus_== 1 ) I_total = S.I_nu[0];

    return I_total;
}

double LOS_data::integrate_LOS_with_binning( int binning_type, int N_bins )
{
    /* Spectral LOS integration using the exact transport eq. sol. for constant property slabs */
    /* ...and using binned spectral coefficients                                                   */

    // 0. Create the spectral bins
    vector<SpectralBin*> B;
    this->create_spectral_bins( binning_type, N_bins, B);

    // 1. Create set of Binned coefficient spectra
    for ( int irp=0; irp<nrps_; irp++ ) {
	rpoints_[irp]->Y_ = new BinnedCoeffSpectra( rpoints_[irp]->X_, B );
    }

    // 2. Create binned intensity spectra for backlight and frontlight
    BinnedSpectralIntensity I_i = BinnedSpectralIntensity( rsm_, T_i_, B );
    BinnedSpectralIntensity I_f = BinnedSpectralIntensity( rsm_, T_f_, B );

    // Spatial integration over discretised points
    for ( int irp=0; irp<nrps_; irp++ ) {
	double ds = rpoints_[irp]->ds_;
	// Spatial integration for each spectral bin
	for ( int iB=0; iB < N_bins; iB++ ) {
	    // exact solution for constant property slab
	    double decay_factor = exp( - ds * rpoints_[irp]->Y_->kappa_bin[iB] );
	    // newly radiated intensity
	    double tmpA = rpoints_[irp]->Y_->j_bin[iB] / rpoints_[irp]->Y_->kappa_bin[iB] * ( 1.0 - decay_factor );
	    if ( isnan(tmpA) ) {
		// tmpA will go NaN when kappa_nu is zero.  Intensity is zero.
		tmpA = 0.0;
	    }
	    // incoming intensity
	    double tmpB = I_i.I_bin[iB] * decay_factor;
	    // write-over old I_bin with new I_bin
	    I_i.I_bin[iB] = tmpA + tmpB;
	}
    }
    // Subtract out Planck intensity from outer BC and cumpute integrated (via summing) intensity
    double I_total = 0.0;
    for ( int iB=0; iB < N_bins; iB++ ) {
	I_i.I_bin[iB] -= I_f.I_bin[iB];
        I_total += I_i.I_bin[iB];
    }

    return I_total;
}

double LOS_data::integrate_LOS_MC( SpectralIntensity &S, int nphotons )
{
    /* Spectral LOS integration using the Monte Carlo method */
    // NOTE: only computing absorption
    // Check vector sizes in S for consistency
    if ( int(S.nu.size()) != nnus_ || int(S.I_nu.size()) != nnus_ ) {
        cout << "Mismatch in size of spectral vector in structure S. " << endl
             << "nnus = " << nnus_ << ", S.nu.size() = " << S.nu.size() << endl
             << "Bailing out!" << endl;
        exit( BAD_INPUT_ERROR );
    }

    double E_total = S.integrate_intensity_spectra();
    double E_out = 0.0;
    
    CRandomMersenne rg = CRandomMersenne((int32)time(0));

    // Loop over all photons
    for ( int ip=0; ip < nphotons; ip++ ) {
        double E_photon = E_total / nphotons;
        double nu = S.random_frequency(rg.Random());
        // int inu = S.random_frequency_interval(rg.Random());
        // cout << "nu = " << nu << endl;
        for ( int irp=0; irp<nrps_; irp++ ) {
            double ds = rpoints_[irp]->ds_;
            double kappa = rpoints_[irp]->X_->kappa_from_nu(nu);
            // double kappa = rpoints_[irp]->X_->kappa_nu[inu];
            double decay_factor = exp( - ds * kappa );
            E_photon *= decay_factor;
        }
        E_out += E_photon;
    }

    return E_out;
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

void LOS_data::write_gas_data_to_file( void )
{
    /* 1. Setup the output file. */
    ofstream ofile;
    ofile.open("LOS-gas-data.txt");
    if( ofile.fail() ) {
	cout << "Error opening file: LOS-gas-data.txt" << endl;
	cout << "Bailing Out!\n";
	exit(FILE_ERROR);
    }
    
    ofile << setprecision(12) << scientific << showpoint;
    
    size_t ntm = rpoints_[0]->Q_->T.size();
    size_t nsp = rpoints_[0]->Q_->massf.size();

    ofile << "# File: LOS-gas-data.txt" << endl
    	  << "# Column 1: Position, s (m)" << endl
          << "# Column 2: Pressure, p (Pa)" << endl
          << "# Column 3: Electron pressure, p_e (Pa)" << endl
          << "# Column 4: Density, rho (kg/m**3" << endl;
    for ( size_t itm=0; itm<ntm; ++itm )
    	ofile << "# Column " << itm + 5 << ": Temperature, T[" << itm << "] (K)" << endl;
    for ( size_t isp=0; isp<nsp; ++isp )
    	ofile << "# Column " << isp + 5 + ntm << ": Mass-fraction, massf[" << isp << "] (ND)" << endl;
           
    for ( int irp=0; irp<nrps_; ++irp ) {
    	RadiatingPoint * rp = this->get_rpoint_pointer(irp);
    	ofile << setw(20) << rp->s_
    	      << setw(20) << rp->Q_->p
    	      << setw(20) << rp->Q_->p_e
    	      << setw(20) << rp->Q_->rho;
    	for ( size_t itm=0; itm<ntm; ++itm )
    	    ofile << setw(20) << rp->Q_->T[itm];
    	for ( size_t isp=0; isp<nsp; ++isp )
    	    ofile << setw(20) << rp->Q_->massf[isp];
    	ofile << endl;
    }
    
    ofile.close();
    	
    return;
}

void LOS_data::load_spectra( int inu_start, int inu_end )
{
    for ( int irp=0; irp<nrps_; ++irp ) {
        RadiatingPoint * rp = this->get_rpoint_pointer(irp);
        ostringstream oss;
        oss << "rpoint_" << irp << ".txt";
        rp->X_->read_from_file( oss.str(), inu_start, inu_end );
    }
}

/* ---------- Tangent-slab data storage and manipulation ----------- */

TS_data::TS_data( RadiationSpectralModel * rsm, int nrps, double T_i, double T_f, bool store_spectra_on_disk )
 : LOS_data( rsm, nrps, T_i, T_f, store_spectra_on_disk )
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

double TS_data::solve_for_divq_OT()
{
    // Get the optically thin solution

    double q_int, q_plus = 0.0;
    
    /* ------- Shock-to-wall (plus) direction ------- */

    // Outer-boundary blackbody intensity (note optical thickness = 0)
    q_int = 0.0;
    for ( int inu=0; inu < nnus_; inu++ ) {
        F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_i_ ) * E_3( 0.0 );
        // Frequency integration
        if ( inu>0 )
            q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
    }

    // Add integrated flux
    q_plus += q_int;

    // get optically thin emission for the slabs
    for ( int irp=0; irp < nrps_; irp++ ) {
        q_int = 2.0 * M_PI * rpoints_[irp]->ds_ * rpoints_[irp]->X_->integrate_emission_spectra();
        *(rpoints_[irp]->Q_rE_rad_) = 4.0 * M_PI * rpoints_[irp]->X_->integrate_emission_spectra();
        // Add integrated flux
        q_plus += q_int;
        // Add contribution to the flux spectra
        for ( int inu=0; inu < nnus_; inu++ )
            F_->q_nu[inu] += 2.0 * M_PI * rpoints_[irp]->ds_ * rpoints_[irp]->X_->j_nu[inu];
    }

    return q_plus;
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
            // The code below can be used to confirm that the E_3() function is returning sensible results
#           if 0
            if ( E_3( tau_nu ) > 0.5 ) {
                cout << "E_3( tau_nu = " << tau_nu << " ) = " << E_3( tau_nu ) << endl;
                exit(0);
            }
#           endif
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

            // The code below can be used to confirm that the E_3() function is returning sensible results
#           if 0
            if ( E_3( tau_nu ) > 0.5 ) {
                cout << "E_3( tau_nu = " << tau_nu << " ) = " << E_3( tau_nu ) << endl;
                exit(0);
            }
#           endif
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

double TS_data::quick_solve_for_divq_in_blocks( int nblks )
{
    // Simple method from Sebastian Karl's VKI thesis where q+ and q- between
    // each node location are solved for
    // This version does the calculation in multiple blocks in an attempt to reduce the memory requirements
    if ( !store_spectra_on_disk_ ) {
        cout << "TS_data::quick_solve_for_divq_in_blocks()" << endl
             << "Currently only implemented for spectra stored on disk" << endl;
        exit(FAILURE);
    }

    double sps_av = double(nnus_) / double(nblks);

    // Make vectors to store integrated fluxes
    vector<double> q_plus(nrps_+1), q_minus(nrps_+1);

    /* Loop over blocks */
    for ( int iblk=0; iblk<nblks; ++iblk ) {
        // Calculate inu_lower and inu_upper
        int nnus_star = 0;
        if ( iblk==0 ) {
            nnus_star = nnus_ - ( nblks - 1 )*int(sps_av);
        }
        else {
            nnus_star = int(sps_av);
        }
        int inu_lower = 0;
        if ( iblk>0 ) inu_lower = nnus_ - ( nblks - 1 )*int(sps_av) + (iblk-1)*int(sps_av);
        int inu_upper = inu_lower + nnus_star - 1;

        cout << "nnus = " << nnus_ << ", nnus_star = " << nnus_star << ", iblk = " << iblk << ", inu_lower = " << inu_lower << ", inu_upper = " << inu_upper << endl;

        /* load the spectra into memory */
        if ( store_spectra_on_disk_ )
            this->load_spectra( inu_lower, inu_upper );

        /* ------- Wall-to-shock (minus) direction ------- */
        int iqm=0;

        // Outer-boundary blackbody intensity (note optical thickness = 0)
        double q_int = 0.0;
        for ( int inu=inu_lower; inu <=inu_upper; inu++ ) {
            F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_f_ ) * E_3( 0.0 );
            // Frequency integration
            if ( inu>0 )
                q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
        }

        // Save integrated flux
        q_minus[iqm] += q_int; iqm++;

        // integrate over emitting slabs using trapezoidal method
        for ( int irp=(nrps_-1); irp >= 0; --irp ) {
            // cout << "Spatial location: | " << irp << " | \r" << flush;
            q_int = 0.0;
            for ( int inu=inu_lower; inu<=inu_upper; inu++ ) {
                // emission and absorption coefficients for this slab
                double kappa_nu = rpoints_[irp]->X_->kappa_nu[inu-inu_lower];
                double j_nu = rpoints_[irp]->X_->j_nu[inu-inu_lower];
                // calculate optical thickness between intervals
                double tau_nu = fabs(rpoints_[irp]->ds_) * kappa_nu;
                // Attenuate incident flux )
                F_->q_nu[inu] *= 2.0 * E_3( tau_nu );
                // Add contribution from this slab
                double tmpA = 0.0;
                if ( j_nu != 0.0 )  tmpA = j_nu / kappa_nu;
                F_->q_nu[inu] += tmpA * M_PI * ( 1.0 -  2.0 * E_3( tau_nu ) );
                // Integrate
                if ( inu>0 )
                q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
            }
            // Save integrated flux
            q_minus[iqm] += q_int; iqm++;
        }

        /* ------- Shock-to-wall (plus) direction ------- */
        int iqp=0;

        // Outer-boundary blackbody intensity (note optical thickness = 0)
        q_int = 0.0;
        for ( int inu=inu_lower; inu<=inu_upper; inu++ ) {
            F_->q_nu[inu] = 2.0 * M_PI * planck_intensity( F_->nu[inu], T_i_ ) * E_3( 0.0 );
            // Frequency integration
            if ( inu>0 )
                q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
        }

        // Save integrated flux
        q_plus[iqp] += q_int; iqp++;

        // integrate over emitting slabs using trapezoidal method
        for ( int irp=0; irp < nrps_; irp++ ) {
            // cout << "Spatial location: | " << irp << " | \r" << flush;
            q_int = 0.0;
            for ( int inu=inu_lower; inu<=inu_upper; inu++ ) {
                // emission and absorption coefficients for this slab
                double kappa_nu = rpoints_[irp]->X_->kappa_nu[inu-inu_lower];
                double j_nu = rpoints_[irp]->X_->j_nu[inu-inu_lower];
                // calculate optical thickness between intervals
                double tau_nu = fabs(rpoints_[irp]->ds_) * kappa_nu;
                // Attenuate incident flux )
                F_->q_nu[inu] *= 2.0 * E_3( tau_nu );
                // Add contribution from this slab
                double tmpA = 0.0;
                if ( j_nu != 0.0 )  tmpA = j_nu / kappa_nu;
                F_->q_nu[inu] += tmpA * M_PI * ( 1.0 -  2.0 * E_3( tau_nu ) );
                // Integrate
                if ( inu>0 )
                    q_int += 0.5 * ( F_->q_nu[inu] + F_->q_nu[inu-1] ) * fabs( F_->nu[inu] - F_->nu[inu-1] );
            }
            // Save integrated flux
            q_plus[iqp] += q_int; iqp++;
        }
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

double TS_data::quick_solve_for_divq_with_binning( int binning_type, int N_bins )
{
    // Simple method from Sebastian Karl's VKI thesis where q+ and q- between
    // each node location are solved for
    // ...with spectral binning

    // Make vectors to store integrated fluxes
    vector<double> q_plus, q_minus;

    // 0. Create the spectral bins
    vector<SpectralBin*> B;
    this->create_spectral_bins( binning_type, N_bins, B);

    // 1. Create set of Binned coefficient spectra
    for ( int irp=0; irp<nrps_; irp++ ) {
	rpoints_[irp]->Y_ = new BinnedCoeffSpectra( rpoints_[irp]->X_, B );
    }

    /* ------- Wall-to-shock (minus) direction ------- */

    // Outer-boundary blackbody intensity (note optical thickness = 0)
    BinnedSpectralFlux F_f( rsm_, T_f_, B );
    double q_int = 0.0;
    for ( int iB=0; iB < N_bins; iB++ ) {
	q_int += F_f.q_bin[iB];
    }

    // Save integrated flux
    q_minus.push_back( q_int );

    // integrate over emitting slabs using trapezoidal method
    for ( int irp=(nrps_-1); irp >= 0; --irp ) {
        // cout << "Spatial location: | " << irp << " | \r" << flush;
        q_int = 0.0;
        for ( int iB=0; iB < N_bins; iB++ ) {
            // emission and absorption coefficients for this slab
            double kappa = rpoints_[irp]->Y_->kappa_bin[iB];
            double j = rpoints_[irp]->Y_->j_bin[iB];
            // calculate optical thickness between intervals
            double tau = fabs(rpoints_[irp]->ds_) * kappa;
            // Attenuate incident flux )
            F_f.q_bin[iB] *= 2.0 * E_3( tau );
            // Add contribution from this slab
	    double tmpA = 0.0;
	    if ( j != 0.0 )	tmpA = j / kappa;
	     F_f.q_bin[iB] += tmpA * M_PI * ( 1.0 -  2.0 * E_3( tau ) );
            // Integrate
            q_int += F_f.q_bin[iB];
        }
        // Save integrated flux
        q_minus.push_back( q_int );
    }

    /* ------- Shock-to-wall (plus) direction ------- */

    // Outer-boundary blackbody intensity (note optical thickness = 0)
    BinnedSpectralFlux F_i( rsm_, T_i_, B );
    q_int = 0.0;
    for ( int iB=0; iB < N_bins; iB++ ) {
	q_int += F_i.q_bin[iB];
    }

    // Save integrated flux
    q_plus.push_back( q_int );

    // integrate over emitting slabs using trapezoidal method
    for ( int irp=0; irp < nrps_; irp++ ) {
        // cout << "Spatial location: | " << irp << " | \r" << flush;
        q_int = 0.0;
        for ( int iB=0; iB < N_bins; iB++ ) {
            // emission and absorption coefficients for this slab
            double kappa = rpoints_[irp]->Y_->kappa_bin[iB];
            double j = rpoints_[irp]->Y_->j_bin[iB];
            // calculate optical thickness between intervals
            double tau = fabs(rpoints_[irp]->ds_) * kappa;
            // Attenuate incident flux )
            F_i.q_bin[iB] *= 2.0 * E_3( tau );
            // Add contribution from this slab
	    double tmpA = 0.0;
	    if ( j != 0.0 )	tmpA = j / kappa;
	    F_i.q_bin[iB] += tmpA * M_PI * ( 1.0 -  2.0 * E_3( tau ) );
            // Integrate
            q_int += F_i.q_bin[iB];
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
		ds_i = rpoints_[irp]->ds_;
		// if ( irp != 0 )
		//    ds_i = 2.0 * ( rpoints_[irp-1]->s_ - rpoints_[irp]->s_ - 0.5 * ds_i );
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

double TS_data::quick_solve_for_planck_divq( double kappa_nu )
{
    // Simple method from Sebastian Karl's VKI thesis where q+ and q- between
    // each node location are solved for
    // using planck function to generate emission coefficients

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
            // calculate optical thickness between intervals
            double tau_nu = fabs(rpoints_[irp]->ds_) * kappa_nu;
            // Attenuate incident flux )
            F_->q_nu[inu] *= 2.0 * E_3( tau_nu );
            // Add contribution from this slab
             F_->q_nu[inu] += 2.0 * M_PI * planck_intensity( F_->nu[inu], rpoints_[irp]->Q_->T[0] ) * E_3( 0.0 );
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
            // calculate optical thickness between intervals
            double tau_nu = fabs(rpoints_[irp]->ds_) * kappa_nu;
            // Attenuate incident flux )
            F_->q_nu[inu] *= 2.0 * E_3( tau_nu );
            // Add contribution from this slab
             F_->q_nu[inu] += 2.0 * M_PI * planck_intensity( F_->nu[inu], rpoints_[irp]->Q_->T[0] ) * E_3( 0.0 );
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
