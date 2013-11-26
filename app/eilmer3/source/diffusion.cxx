/** \file cns_diffusion.cxx
 *  \ingroup mbcns2
 *  \brief Definitions for the diffusion model class and functions.
 **/

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "kernel.hh"
#include "diffusion.hh"

using namespace std;

DiffusionModel::DiffusionModel(const string name, int nsp)
    : name_(name), nsp_(nsp), e_index_( -1 )
{
    Gas_model *gmodel = get_gas_model_ptr();
    min_molef_ = DEFAULT_MIN_MOLE_FRACTION;	// use default from gas-models
    small_DAB_ = 1.0e-20;
    DAV_im_.resize( nsp_ );
    x_.resize( nsp_ );
    M_.resize( nsp_ );
    Z_.resize( nsp_ );
    for( int isp = 0; isp < nsp_; ++isp ) {
	M_[isp] = gmodel->molecular_weight(isp);
	if ( gmodel->species_name(isp)=="e_minus" ) e_index_ = isp;
	if ( gmodel->species_name(isp).find("minus")!=string::npos ) Z_[isp] = -1.0;
	else if ( gmodel->species_name(isp).find("plus")!=string::npos ) Z_[isp] = 1.0;
	else Z_[isp] = 0.0;
    }
}

DiffusionModel::DiffusionModel(const DiffusionModel &d)
    : name_(d.name_), nsp_(d.nsp_), e_index_( d.e_index_ ), min_molef_( d.min_molef_ ),
      small_DAB_( d.small_DAB_ ), DAV_im_( d.DAV_im_ ), x_( d.x_ ), M_( d.M_ ),
      Z_( d.Z_ )
{}

DiffusionModel::~DiffusionModel() {}

void
DiffusionModel::
fill_in_x(double rho, const vector<double> &massf)
{
    // Compute total moles per kg
    double total_moles_per_kg = 0.0;
    for ( size_t isp = 0; isp < massf.size(); ++isp )
    	total_moles_per_kg += massf[isp]/M_[isp];
    
    // Compute mole fractions
    for ( size_t isp = 0; isp < massf.size(); ++isp )
	x_[isp] = massf[isp]/M_[isp]/total_moles_per_kg;
}

void
DiffusionModel::
fill_in_DAV_im(const matrix &D_AB)
{
    // 1. Heavy particles
    for ( int isp = 0; isp < nsp_; ++isp ) {
#       if WILKE_MIXING_RULE_WITH_AMBIPOLAR_CORRECTION
    	if ( isp==e_index_ ) continue;
#       endif
	double sum = 0.0;
	for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	    if ( isp == jsp ) continue;
	    // The following two if-statements should generally catch the
	    // same flow condition, namely, a zero or very small presence of
	    // a certain species.  In this case the diffusion is effectively
	    // zero and its contribution to the mixture diffusion coefficient
	    // may be ignored.
	    //
	    // The two statements are used for extra security in detecting the
	    // condition.
	    if ( D_AB[isp][jsp] <= small_DAB_ ) continue;  // there is effectively nothing to diffuse
	    if ( x_[jsp] > min_molef_ ) {
		sum += x_[jsp] / D_AB[isp][jsp];
	    }
	}
	if ( sum <= 0.0 ) {
	    DAV_im_[isp] = 0.0;
	}
	else {
	    DAV_im_[isp] = (1.0 - x_[isp]) / sum;
	}
    }
    
#   if WILKE_MIXING_RULE_WITH_AMBIPOLAR_CORRECTION
    // 2. Free electrons and ions with ambipolar correction
    if ( e_index_ > -1 ) {
    	double Dax_ion_sum = 0.0;
    	double Mx_ion_sum = 0.0;
    	for ( int isp = 0; isp < nsp_; ++isp ) {
    	    if ( Z_[isp] > 0.0 ) {
    	    	DAV_im_[isp] *= 2.0;
    	    	Dax_ion_sum += DAV_im_[isp] * x_[isp];
    	    	Mx_ion_sum += M_[isp] * x_[isp];
    	    }
    	}
    	DAV_im_[e_index_] = M_[e_index_] * Dax_ion_sum / Mx_ion_sum;
    	if ( !isfinite( DAV_im_[e_index_] ) ) {
    	    // Dax_ion_sum and Mx_ion_sum were probably zero
    	    DAV_im_[e_index_] = 0.0;
    	}
    }
#   endif
}

NoDiffusionModel::
NoDiffusionModel(const string name, int nsp)
    : DiffusionModel(name, nsp) {}

NoDiffusionModel::
NoDiffusionModel(const NoDiffusionModel &s)
    : DiffusionModel(s.name_, s.nsp_)
{}

NoDiffusionModel::
~NoDiffusionModel() {}

string
NoDiffusionModel::
str() const
{
    ostringstream ost;
    cout << "NoDiffusionModel::str() - not implemented.\n";
    ost << "NoDiffusionModel";
    return ost.str();
}


void
NoDiffusionModel::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const std::vector<double> &dfdx, 
			   const std::vector<double> &dfdy,
			   const std::vector<double> &dfdz,
			   std::vector<double> &jx, 
			   std::vector<double> &jy,
			   std::vector<double> &jz)
{
    // Just ensure all fluxes are 0.
    for( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = 0.0;
	jy[isp] = 0.0;
	jz[isp] = 0.0;
    }
}

FicksFirstLaw::
FicksFirstLaw(const string name, int nsp)
    : DiffusionModel(name, nsp)
{}

FicksFirstLaw::
FicksFirstLaw(const FicksFirstLaw &f)
    : DiffusionModel(f.name_, f.nsp_)
{}

FicksFirstLaw::
~FicksFirstLaw() {}

string
FicksFirstLaw::
str() const
{
    ostringstream ost;
    cout << "FicksFirstLaw::str() - not implemented.\n";
    ost << "FicksFirstLaw";
    return ost.str();
}

void
FicksFirstLaw::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    fill_in_x(Q.rho, Q.massf);
    fill_in_DAV_im(Q.D_AB);

    // Set diffusive fluxes...
    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
	jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
	jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }

#   if FICKS_WITH_SUTTON_AND_GNOFFO_CORRECTION
    // Correction as suggested by Sutton and Gnoffo, 1998  
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;
    
    for ( int isp = 0; isp < nsp_; ++isp ) {
	sum_x += jx[isp];
	sum_y += jy[isp];
	sum_z += jz[isp];
    }

    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = jx[isp] - Q.massf[isp] * sum_x;
	jy[isp] = jy[isp] - Q.massf[isp] * sum_y;
	jz[isp] = jz[isp] - Q.massf[isp] * sum_z;
    }
#   endif

}

StefanMaxwellModel::
StefanMaxwellModel(const string name, int nsp, int no_iters)
    : DiffusionModel(name, nsp), no_iters_(no_iters)
{
    j0_.resize( nsp_ );
    j1_.resize( nsp_ );
}

StefanMaxwellModel::
StefanMaxwellModel( const StefanMaxwellModel &s )
    : DiffusionModel( s.name_, s.nsp_ ), no_iters_( s.no_iters_ )
{
    j0_.resize( nsp_ );
    j1_.resize( nsp_ );
}

StefanMaxwellModel::
~StefanMaxwellModel() {}

string
StefanMaxwellModel::
str() const
{
    ostringstream ost;
    cout << "StefanMaxwellModel::str() - not implemented.\n";
    ost << "StefanMaxwellModel";
    return ost.str();
}

void
StefanMaxwellModel::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    // NOTE: Presently the Stefan-Maxwell model is not compatible
    // with modelled turbulence.

    UNUSED_VARIABLE(D_t);

    double sum = 0.0;
    const double tol = 1.0e-6; // See Sutton And Gnoffo (1998) p.9
    bool converged = true;

    fill_in_x(Q.rho, Q.massf);
    fill_in_DAV_im(Q.D_AB);

    // Inline calculation of mixture molecular weight
    double temp = 0.0;
    for ( int isp = 0; isp < nsp_; ++isp ) {
	temp += Q.massf[isp]/M_[isp];
    }
    double MW_mix = 1.0/temp;

    // -----------------------------------
    // 1. In the x-direction...
    // -----------------------------------

    // Establish a guess in j0_;
    for( int isp = 0; isp < nsp_; ++isp ) {
	j0_[isp] = -Q.rho * DAV_im_[isp] * dfdx[isp];
    }

    for( int iter = 0; iter < no_iters_; ++iter ) {
	converged = true;
	for( int isp = 0; isp < nsp_; ++isp ) {
	    sum = 0.0;
	    if( DAV_im_[isp] <= 0.0 ) {
		// No diffusion because it is completely this species.
		j1_[isp] = 0.0;
		continue;
	    }
	    for( int jsp = 0; jsp < nsp_; ++jsp ) {
		if( isp == jsp) continue;
		
		sum += Q.rho * (MW_mix / M_[jsp]) * dfdx[jsp] +
		    (MW_mix / M_[jsp]) * (j0_[jsp] / Q.D_AB[isp][jsp]);
		
	    }
	    j1_[isp] = -Q.rho * DAV_im_[isp] * dfdx[isp] + 
		Q.massf[isp] * DAV_im_[isp] / (1.0 - x_[isp] ) * sum;
	}
	sum = 0.0;
	for( int jsp = 0; jsp < nsp_; ++jsp ) {
	    sum += j1_[jsp];
	}
	for( int isp = 0; isp < nsp_; ++isp ) {
	    jx[isp] = j1_[isp] - Q.massf[isp] * sum;
	    if( fabs(j1_[isp] - jx[isp]) > tol ) {
		converged = false;
	    }
	    j0_[isp] = jx[isp];
	    
	}
	if( converged ) {
	    //cout << "no iterations for jx= " << iter << endl;
	    break;
	}
    }
    
    // -----------------------------------------------
    // 2. in the y-direction
    // -----------------------------------------------

    // Establish a guess in j0_;
    for( int isp = 0; isp < nsp_; ++isp ) {
	j0_[isp] = -Q.rho * DAV_im_[isp] * dfdy[isp];
    }

    for( int iter = 0; iter < no_iters_; ++iter ) {
	converged = true;
	for( int isp = 0; isp < nsp_; ++isp ) {
	    sum = 0.0;
	    if( DAV_im_[isp] <= 0.0 ) {
		// No diffusion because it is completely this species.
		j1_[isp] = 0.0;
		continue;
	    }
	    for( int jsp = 0; jsp < nsp_; ++jsp ) {
		if( isp == jsp) continue;
		
		sum += Q.rho * (MW_mix / M_[jsp]) * dfdy[jsp] +
		    (MW_mix / M_[jsp]) * (j0_[jsp] / Q.D_AB[isp][jsp]);
	    }
	    j1_[isp] = -Q.rho * DAV_im_[isp] * dfdy[isp] + 
		Q.massf[isp] * DAV_im_[isp] / (1.0 - x_[isp] ) * sum;
	}
	sum = 0.0;
	for( int jsp = 0; jsp < nsp_; ++jsp ) {
	    sum += j1_[jsp];
	}
	for( int isp = 0; isp < nsp_; ++isp ) {
	    jy[isp] = j1_[isp] - Q.massf[isp] * sum;
	    if( fabs(j1_[isp] - jy[isp]) > tol ) {
		converged = false;
	    }
	    j0_[isp] = jy[isp];
	}
	if( converged ) {
	    break;
	}
    }

    // -----------------------------------------------
    // 3. in the z-direction
    // -----------------------------------------------

    // Establish a guess in j0_;
    for( int isp = 0; isp < nsp_; ++isp ) {
	j0_[isp] = -Q.rho * DAV_im_[isp] * dfdz[isp];
    }

    for( int iter = 0; iter < no_iters_; ++iter ) {
	converged = true;
	for( int isp = 0; isp < nsp_; ++isp ) {
	    sum = 0.0;
	    if( DAV_im_[isp] <= 0.0 ) {
		// No diffusion because it is completely this species.
		j1_[isp] = 0.0;
		continue;
	    }
	    for( int jsp = 0; jsp < nsp_; ++jsp ) {
		if( isp == jsp) continue;
		
		sum += Q.rho * (MW_mix / M_[jsp]) * dfdz[jsp] +
		    (MW_mix / M_[jsp]) * (j0_[jsp] / Q.D_AB[isp][jsp]);
	    }
	    j1_[isp] = -Q.rho * DAV_im_[isp] * dfdz[isp] + 
		Q.massf[isp] * DAV_im_[isp] / (1.0 - x_[isp] ) * sum;
	}
	sum = 0.0;
	for( int jsp = 0; jsp < nsp_; ++jsp ) {
	    sum += j1_[jsp];
	}
	for( int isp = 0; isp < nsp_; ++isp ) {
	    jz[isp] = j1_[isp] - Q.massf[isp] * sum;
	    if( fabs(j1_[isp] - jz[isp]) > tol ) {
		converged = false;
	    }
	    j0_[isp] = jz[isp];
	}
	if( converged ) {
	    break;
	}
    }
}

// Self Consistent Effective Binary Diffusion Method (SCEBDM) of Ramshaw and Chang [JNET Vol 15 1990]

RamshawChangModel::
RamshawChangModel( const string name, int nsp )
    : DiffusionModel( name, nsp )
{
    dcdx_.resize( nsp_ );
    dcdy_.resize( nsp_ );
}

RamshawChangModel::
RamshawChangModel( const RamshawChangModel &rcm )
    : DiffusionModel( rcm.name_, rcm.nsp_ )
{
    dcdx_.resize( nsp_ );
    dcdy_.resize( nsp_ );
}

RamshawChangModel::
~RamshawChangModel()
{}

string
RamshawChangModel::
str() const
{
    ostringstream ost;
    cout << "RamshawChangModel::str() - not implemented.\n";
    ost << "RamshawChangModel";
    return ost.str();
}

void
RamshawChangModel::
fill_in_dcdx_dcdy( double rho, const vector<double> &dfdx,
                               const vector<double> &dfdy )
{
    for ( size_t isp=0; isp<M_.size(); ++isp ) {
	dcdx_[isp] = rho*dfdx[isp]/M_[isp];
	dcdy_[isp] = rho*dfdy[isp]/M_[isp];
    }
    
    return;
}

void
RamshawChangModel::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    fill_in_x(Q.rho, Q.massf);
    fill_in_DAV_im(Q.D_AB);
    fill_in_dcdx_dcdy( Q.rho, dfdx, dfdy );

    //    cout << "rho= " << rho << endl;
    
    // 0. Pre-calculate jsp summation terms
    double MDdcdx = 0.0;
    double MDdcdy = 0.0;
    double MqrhoD = 0.0;
    for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	if ( jsp==e_index_ ) continue;
	MDdcdx += M_[jsp] * DAV_im_[jsp] * dcdx_[jsp];
	MDdcdy += M_[jsp] * DAV_im_[jsp] * dcdy_[jsp];
	MqrhoD += M_[jsp] * Z_[jsp] * Q.rho * Q.massf[jsp] * DAV_im_[jsp];
    }
    
    // A. Heavy species (non-electrons)
    
    // 1.  In the x-direction...
    for( int isp = 0; isp < nsp_; ++isp ) {
        if ( isp==e_index_ || x_[isp] < min_molef_ ) continue;
        double As = 0.0;
        if ( x_[e_index_] >= min_molef_ ) {
	    As = ( 1.0 / Z_[e_index_] * Q.rho * Q.massf[e_index_] ) * 
		 ( M_[isp] * Z_[isp] * Q.rho * Q.massf[isp] * DAV_im_[isp] - Q.massf[isp] * MqrhoD ) * dcdx_[e_index_];
        }
	jx[isp] = -M_[isp] * DAV_im_[isp] * dcdx_[isp] + Q.massf[isp] * MDdcdx + As;
    }
    
    // 2. in the y-direction
    for( int isp = 0; isp < nsp_; ++isp ) {
	if ( isp==e_index_ || x_[isp] < min_molef_ ) continue;
        double As = 0.0;
        if ( x_[e_index_] >= min_molef_ ) {
	    As = ( 1.0 / Z_[e_index_] ) / ( Q.rho * Q.massf[e_index_] ) * 
		 ( M_[isp] * Z_[isp] * Q.rho * Q.massf[isp] * DAV_im_[isp] - Q.massf[isp] * MqrhoD ) * dcdy_[e_index_];
        }
	jy[isp] = - M_[isp] * DAV_im_[isp] * dcdy_[isp] + Q.massf[isp] * MDdcdy + As;
    }
    
    // B. Electrons
    jx[e_index_] = 0.0;
    jy[e_index_] = 0.0;
    if ( x_[e_index_] >= min_molef_ ) {
	for( int isp = 0; isp < nsp_; ++isp ) {
	    if ( isp==e_index_ ) continue;
	    jx[e_index_] += Z_[isp] * jx[isp];
	    jy[e_index_] += Z_[isp] * jy[isp];
	}
	jx[e_index_] /= - Z_[e_index_];
	jy[e_index_] /= - Z_[e_index_];
    }

    // 9. Thermal nonequilibrium (T!=Te) correction
    //    Correction factor = p/RT / c = ( c_HP + c_e * T_e / T ) / c
    //    (see Eq. 22 in Ramshaw et al PCPP Vol. 13 | No. 3 1993)
    if ( Q.T.size()>1 ) {
        double T = Q.T.front();
        double c = 0.0;
        for( int isp = 0; isp < nsp_; ++isp )
            c += x_[isp];
        double f_neq = Q.p/PC_R_u/T/c;
        for( int isp = 0; isp < nsp_; ++isp ) {
            jx[isp] *= f_neq;
            jy[isp] *= f_neq;
        }
    }
    
    return;
}

ConstantLewisNumber::
ConstantLewisNumber(const string name, int nsp)
    : DiffusionModel(name, nsp)
{
    // Just set the Lewis number here for the number
    Le_ = 1.4;
}

ConstantLewisNumber::
ConstantLewisNumber(const ConstantLewisNumber &c)
    : DiffusionModel(c.name_, c.nsp_), Le_( c.Le_ )
{}

ConstantLewisNumber::
~ConstantLewisNumber() {}

string
ConstantLewisNumber::
str() const
{
    ostringstream ost;
    cout << "ConstantLewisNumber::str() - not implemented.\n";
    ost << "ConstantLewisNumber";
    return ost.str();
}

void
ConstantLewisNumber::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    // Calculate the Prandtl number
    Gas_model *gmodel = get_gas_model_ptr();
    double k_total = 0.0;
    for ( size_t itm=0; itm<Q.k.size(); ++itm ) {
    	k_total += Q.k[itm];
    }
    int status;
    double Prandtl = Q.mu * gmodel->Cp(Q, status) / k_total;
    
    // Calculate the species average diffusion coefficients
    for ( int isp = 0; isp < nsp_; ++isp )
        DAV_im_[isp] = Q.mu / Q.rho / Prandtl / Le_;

    // NOTE: previously this relation between D and Le was erroneously written as:
    // DAV_im_[isp] = Le_ * Q.mu / Prandtl;
    
#   if 0
    // Apply ambipolar diffusion correction
    if ( e_index_ > -1 ) {
    	double Dax_ion_sum = 0.0;
    	double Mx_ion_sum = 0.0;
    	for ( int isp = 0; isp < nsp_; ++isp ) {
    	    if ( Z_[isp] > 0.0 ) {
    	    	DAV_im_[isp] *= 2.0;
    	    	Dax_ion_sum += DAV_im_[isp] * x_[isp];
    	    	Mx_ion_sum += M_[isp] * x_[isp];
    	    }
    	}
    	DAV_im_[e_index_] = M_[e_index_] * Dax_ion_sum / Mx_ion_sum;
    	if ( !isfinite( DAV_im_[e_index_] ) ) {
    	    // Dax_ion_sum and Mx_ion_sum were probably zero
    	    DAV_im_[e_index_] = 0.0;
    	}
    }
#   endif
    
    // Set diffusive fluxes via Fick's first law
    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
	jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
	jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }
    

    	
    
    return;
}

ConstantSchmidtNumber::
ConstantSchmidtNumber(const string name, int nsp)
    : DiffusionModel(name, nsp)
{
    // Just set the Schmidt number here for the moment
    Sc_ = 0.7;
}

ConstantSchmidtNumber::
ConstantSchmidtNumber(const ConstantSchmidtNumber &c)
    : DiffusionModel(c.name_, c.nsp_), Sc_( c.Sc_ )
{}

ConstantSchmidtNumber::
~ConstantSchmidtNumber() {}

string
ConstantSchmidtNumber::
str() const
{
    ostringstream ost;
    cout << "ConstantSchmidtNumber::str() - not implemented.\n";
    ost << "ConstantSchmidtNumber";
    return ost.str();
}

void
ConstantSchmidtNumber::
calculate_diffusion_fluxes(const Gas_data &Q,
                           double D_t,
                           const vector<double> &dfdx,
                           const vector<double> &dfdy,
                           const vector<double> &dfdz,
                           vector<double> &jx,
                           vector<double> &jy,
                           vector<double> &jz)
{
    // Calculate the species average diffusion coefficients
    for ( int isp = 0; isp < nsp_; ++isp )
        DAV_im_[isp] = Q.mu / Q.rho / Sc_;

#   if 0
    // Apply ambipolar diffusion correction
    if ( e_index_ > -1 ) {
        double Dax_ion_sum = 0.0;
        double Mx_ion_sum = 0.0;
        for ( int isp = 0; isp < nsp_; ++isp ) {
            if ( Z_[isp] > 0.0 ) {
                DAV_im_[isp] *= 2.0;
                Dax_ion_sum += DAV_im_[isp] * x_[isp];
                Mx_ion_sum += M_[isp] * x_[isp];
            }
        }
        DAV_im_[e_index_] = M_[e_index_] * Dax_ion_sum / Mx_ion_sum;
        if ( !isfinite( DAV_im_[e_index_] ) ) {
            // Dax_ion_sum and Mx_ion_sum were probably zero
            DAV_im_[e_index_] = 0.0;
        }
    }
#   endif

    // Set diffusive fluxes via Fick's first law
    for ( int isp = 0; isp < nsp_; ++isp ) {
        jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
        jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
        jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }




    return;
}

static DiffusionModel* dmodel = 0;

int set_diffusion_model( const string diffusion_model )
{
    delete dmodel; // In case we're already pointing.
    int nsp = get_gas_model_ptr()->get_number_of_species();

    if( diffusion_model == "Stefan-Maxwell" ) {
	dmodel = new StefanMaxwellModel("Stefan-Maxwell diffusion", nsp, 10 );
    }
    else if ( diffusion_model == "FicksFirstLaw" ) {
	dmodel = new FicksFirstLaw("Fick's first law (of diffusion)", nsp );
    }
    else if ( diffusion_model == "None" ) {
	dmodel = new NoDiffusionModel("No Diffusion", nsp );
    }
    else if ( diffusion_model == "Ramshaw-Chang" ) {
	dmodel = new RamshawChangModel("Ramshaw-Chang diffusion", nsp );
    }
    else if ( diffusion_model == "ConstantLewisNumber" ) {
	dmodel = new ConstantLewisNumber("Constant Lewis number diffusion", nsp );
    }
    else if ( diffusion_model == "ConstantSchmidtNumber" ) {
        dmodel = new ConstantSchmidtNumber("Constant Schmidt number diffusion", nsp );
    }
    else {
	cout << "set_diffusion_model(): " << diffusion_model
	     << " is unknown, bailing out!\n";
	exit( BAD_INPUT_ERROR );
    }
    return 0;

}

void calculate_diffusion_fluxes(const Gas_data &Q,
				double D_t,
				const vector<double> &dfdx, 
				const vector<double> &dfdy,
				const vector<double> &dfdz,
				vector<double> &jx, 
				vector<double> &jy,
				vector<double> &jz)
{
    dmodel->calculate_diffusion_fluxes(Q, D_t, dfdx, dfdy, dfdz, jx, jy, jz);
}

