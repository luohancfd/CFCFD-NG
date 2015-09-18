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

constexpr bool WILKE_MIXING_RULE_WITH_AMBIPOLAR_CORRECTION = true;
constexpr bool FICKS_WITH_SUTTON_AND_GNOFFO_CORRECTION = true;

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
	if ( WILKE_MIXING_RULE_WITH_AMBIPOLAR_CORRECTION ) {
	    if ( isp==e_index_ ) continue;
	}
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
    
    if ( WILKE_MIXING_RULE_WITH_AMBIPOLAR_CORRECTION ) {
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
    }
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

    if ( FICKS_WITH_SUTTON_AND_GNOFFO_CORRECTION ) {
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
    }

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
{}

RamshawChangModel::
RamshawChangModel( const RamshawChangModel &rcm )
    : DiffusionModel( rcm.name_, rcm.nsp_ )
{}

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

    // -1. Calculate DAV_im (effective binary diffusivities)
    // Ramshaw and Chang, Plasma Chem. Plasma Process. 13, 489 (1993) equation (23)
    double summ = 0.0;
    for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	if ( jsp == e_index_ ) continue;
	summ += Q.massf[jsp] / sqrt(M_[jsp]);
    }
    for ( int isp = 0; isp < nsp_; ++isp ) {
	if ( isp == e_index_ ) continue;
	double sum = 0.0;
	for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	    if ( jsp == isp ) continue;
	    if ( jsp == e_index_ ) continue;
	    if ( Q.D_AB[isp][jsp] <= small_DAB_ or isnan(Q.D_AB[isp][jsp]) ) continue;
	    // there is effectively nothing to diffuse
	    sum += x_[jsp] / Q.D_AB[isp][jsp];
	}
	if ( sum <= 0.0 or summ <= 0.0 ) {
	    DAV_im_[isp] = 0.0;
	}
	else {
	    //DAV_im_[isp] = (1.0 - x_[isp]) / sum;
	    DAV_im_[isp] = (1.0 - Q.massf[isp] / sqrt(M_[isp]) / summ) / sum;
	}
    }
    
    // 0. Pre-calculate jsp summation terms
    double MDdcdx = 0.0;
    double MDdcdy = 0.0;
    double MDdcdz = 0.0;
    double MqrhoD = 0.0;
    for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	if ( jsp==e_index_ ) continue;
	MDdcdx += Q.rho * DAV_im_[jsp] * dfdx[jsp];
	MDdcdy += Q.rho * DAV_im_[jsp] * dfdy[jsp];
	MDdcdz += Q.rho * DAV_im_[jsp] * dfdz[jsp];
	MqrhoD += Z_[jsp] * Q.rho * Q.massf[jsp] * DAV_im_[jsp];
    }
    
    // A. Heavy species (non-electrons)

    for( int isp = 0; isp < nsp_; ++isp ) {
        if ( isp==e_index_ ) continue;
        double As = 0.0;
//and G.electric_field_work 
	if ( Q.massf[e_index_]>0.0 and x_[e_index_]>=min_molef_ ) {
            As = ( 1.0 / Z_[e_index_] ) / (Q.rho * Q.massf[e_index_] ) * 
		( Z_[isp] * Q.rho * Q.massf[isp] * DAV_im_[isp] - Q.massf[isp] * MqrhoD ) * 
		Q.rho * dfdx[e_index_];
	}
	jx[isp] = -Q.rho * DAV_im_[isp] * dfdx[isp] + Q.massf[isp] * MDdcdx + As;
        As = 0.0;
	if ( Q.massf[e_index_]>0.0 and x_[e_index_]>=min_molef_ ) {
            As = ( 1.0 / Z_[e_index_] ) / (Q.rho * Q.massf[e_index_] ) * 
		( Z_[isp] * Q.rho * Q.massf[isp] * DAV_im_[isp] - Q.massf[isp] * MqrhoD ) * 
		Q.rho * dfdy[e_index_];
	}
	jy[isp] = -Q.rho * DAV_im_[isp] * dfdy[isp] + Q.massf[isp] * MDdcdy + As;
        As = 0.0;
	if ( Q.massf[e_index_]>0.0 and x_[e_index_]>=min_molef_ ) {
            As = ( 1.0 / Z_[e_index_] ) / (Q.rho * Q.massf[e_index_] ) * 
		( Z_[isp] * Q.rho * Q.massf[isp] * DAV_im_[isp] - Q.massf[isp] * MqrhoD ) * 
		Q.rho * dfdz[e_index_];
	}
	jz[isp] = -Q.rho * DAV_im_[isp] * dfdz[isp] + Q.massf[isp] * MDdcdz + As;
    }
    
    
    // B. Electrons
    if (e_index_ >= 0) {
	jx[e_index_] = 0.0;
	jy[e_index_] = 0.0;
	jz[e_index_] = 0.0;
	for( int isp = 0; isp < nsp_; ++isp ) {
	    if ( isp==e_index_ ) continue;
	    jx[e_index_] += Z_[isp] * jx[isp] / M_[isp];
	    jy[e_index_] += Z_[isp] * jy[isp] / M_[isp];
	    jz[e_index_] += Z_[isp] * jz[isp] / M_[isp];
	}
	jx[e_index_] *= - M_[e_index_] / Z_[e_index_];
	jy[e_index_] *= - M_[e_index_] / Z_[e_index_];
	jz[e_index_] *= - M_[e_index_] / Z_[e_index_];
    }

    // 9. Thermal nonequilibrium (T!=Te) correction
    //    Correction factor = p/RT / c = ( c_HP + c_e * T_e / T ) / c
    //    (see Eq. 22 in Ramshaw et al PCPP Vol. 13 | No. 3 1993)
    if ( Q.T.size()>1 ) {
        double T = Q.T.front();
        double c = 0.0;
        for( int isp = 0; isp < nsp_; ++isp )
            c += Q.rho * Q.massf[isp] / M_[isp];
        double f_neq = Q.p/PC_R_u/T/c;
        for( int isp = 0; isp < nsp_; ++isp ) {
            jx[isp] *= f_neq;
            jy[isp] *= f_neq;
            jz[isp] *= f_neq;
        }
    }
    
    return;
}

ConstantLewisNumber::
ConstantLewisNumber(const string name, int nsp, double Le)
    : DiffusionModel(name, nsp), Le_(Le) {}

ConstantLewisNumber::
ConstantLewisNumber(const ConstantLewisNumber &c)
    : DiffusionModel(c.name_, c.nsp_), Le_(c.Le_) {}

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
ConstantSchmidtNumber(const string name, int nsp, double Sc)
    : DiffusionModel(name, nsp), Sc_(Sc) {}

ConstantSchmidtNumber::
ConstantSchmidtNumber(const ConstantSchmidtNumber &c)
    : DiffusionModel(c.name_, c.nsp_), Sc_(c.Sc_) {}

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


ConstantLewisNumber_DRM19::
ConstantLewisNumber_DRM19(const string name, int nsp)
    : DiffusionModel(name, nsp) {}

ConstantLewisNumber_DRM19::
ConstantLewisNumber_DRM19(const ConstantLewisNumber_DRM19 &c)
    : DiffusionModel(c.name_, c.nsp_) {}

ConstantLewisNumber_DRM19::
~ConstantLewisNumber_DRM19() {}

string
ConstantLewisNumber_DRM19::
str() const
{
    ostringstream ost;
    cout << "ConstantLewisNumber_DRM19::str() - not implemented.\n";
    ost << "ConstantLewisNumber_DRM19";
    return ost.str();
}

//Constant Lewis numbers but different for each species, used for DRM19
//reaction mechanism with differential diffusion
void
ConstantLewisNumber_DRM19::
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
    
    //specify individual Lewis number for each species
    //"[0]-H2" "[1]-H" "[2]-O" "[3]-O2" "[4]-OH" "[5]-H2O"
    //"[6]-HO2" "[7]-CH2" "[8]-CH2_S" "[9]-CH3" "[10]-CH4"
    //"[11]-CO" "[12]-CO2" "[13]-CHO" "[14]-CH2O" "[15]-CH3O"
    //"[16]-C2H4" "[17]-C2H5" "[18]-C2H6" "[19]-N2" "[20]-Ar"
    //"[21]-He"
    vector<double> Le_DRM19;
    Le_DRM19.resize(nsp_, 1.0);
    Le_DRM19[0] = 0.317; // H2
    Le_DRM19[1] = 0.189; // H
    Le_DRM19[2] = 0.712; // O
    Le_DRM19[3] = 1.086; // O2
    Le_DRM19[4] = 0.736; // OH
    Le_DRM19[5] = 0.854; // H2O
    Le_DRM19[6] = 1.079; // HO2
    Le_DRM19[7] = 1.023; // CH2
    Le_DRM19[8] = 1.022; // CH2_S
    Le_DRM19[9] = 1.049; // CH3
    Le_DRM19[10] = 1.043; // CH4
    Le_DRM19[11] = 1.171; // CO
    Le_DRM19[12] = 1.404; // CO2
    Le_DRM19[13] = 1.314; // CHO
    Le_DRM19[14] = 1.329; // CH2O
    Le_DRM19[15] = 1.360; // CH3O
    Le_DRM19[16] = 1.402; // C2H4
    Le_DRM19[17] = 1.551; // C2H5
    Le_DRM19[18] = 1.546; // C2H6
    Le_DRM19[19] = 1.152; // N2
    Le_DRM19[20] = 1.173; // Ar
    // Verhoeven et al. Combustion and Flame, 2012
     
    // Calculate the species average diffusion coefficients
    for ( int isp = 0; isp < nsp_; ++isp )
        DAV_im_[isp] = Q.mu / Q.rho / Prandtl / Le_DRM19[isp];
    
    // Set diffusive fluxes via Fick's first law
    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
	jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
	jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }
    
    return;
}


FicksFirstLaw_DRM19::
FicksFirstLaw_DRM19(const string name, int nsp)
    : DiffusionModel(name, nsp)
{}

FicksFirstLaw_DRM19::
FicksFirstLaw_DRM19(const FicksFirstLaw_DRM19 &f)
    : DiffusionModel(f.name_, f.nsp_)
{}

FicksFirstLaw_DRM19::
~FicksFirstLaw_DRM19() {}

string
FicksFirstLaw_DRM19::
str() const
{
    ostringstream ost;
    cout << "FicksFirstLaw_DRM19::str() - not implemented.\n";
    ost << "FicksFirstLaw_DRM19";
    return ost.str();
}

void
FicksFirstLaw_DRM19::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    //first calculate binary diffusion coefficient D_AB
    //according to the Chapman-Enskog relation
    //refer to "The Properties of Gases and Liquids", Reid et al.

    vector<double> d_c;
    d_c.resize(nsp_, 0.0);
    d_c[0] = 2.920; // H2
    d_c[1] = 2.050; // H
    d_c[2] = 2.750; // O
    d_c[3] = 3.458; // O2
    d_c[4] = 2.750; // OH
    d_c[5] = 2.605; // H2O
    d_c[6] = 3.458; // HO2
    d_c[7] = 3.800; // CH2
    d_c[8] = 3.800; // CH2_S
    d_c[9] = 3.800; // CH3
    d_c[10] = 3.746; // CH4
    d_c[11] = 3.650; // CO
    d_c[12] = 3.763; // CO2
    d_c[13] = 3.590; // CHO
    d_c[14] = 3.590; // CH2O
    d_c[15] = 3.690; // CH3O
    d_c[16] = 3.971; // C2H4
    d_c[17] = 4.302; // C2H5
    d_c[18] = 4.302; // C2H6
    d_c[19] = 3.621; // N2
    d_c[20] = 3.330; // Ar
    d_c[21] = 2.576; // He
    //Lennard-Jones collision diameter "sigma" in Angstroms
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    vector<double> e_LJ;
    e_LJ.resize(nsp_, 0.0);
    e_LJ[0] = 38.000; // H2
    e_LJ[1] = 145.000; // H
    e_LJ[2] = 80.000; // O
    e_LJ[3] = 107.400; // O2
    e_LJ[4] = 80.000; // OH
    e_LJ[5] = 572.400; // H2O
    e_LJ[6] = 107.400; // HO2
    e_LJ[7] = 144.000; // CH2
    e_LJ[8] = 144.000; // CH2_S
    e_LJ[9] = 144.000; // CH3
    e_LJ[10] = 141.400; // CH4
    e_LJ[11] = 98.100; // CO
    e_LJ[12] = 244.000; // CO2
    e_LJ[13] = 498.000; // CHO
    e_LJ[14] = 498.000; // CH2O
    e_LJ[15] = 417.000; // CH3O
    e_LJ[16] = 280.800; // C2H4
    e_LJ[17] = 252.300; // C2H5
    e_LJ[18] = 252.300; // C2H6
    e_LJ[19] = 97.530; // N2
    e_LJ[20] = 136.500; // Ar
    e_LJ[21] = 10.200; // He
    //Lennard-Jones potential well depth "epsilon/kb" in Kelvins
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    double Coef_A = 1.06036;
    double Coef_B = 0.15610;
    double Coef_C = 0.19300;
    double Coef_D = 0.47635;
    double Coef_E = 1.03587;
    double Coef_F = 1.52996;
    double Coef_G = 1.76474;
    double Coef_H = 3.89411;
    //Coefficients needed for calculating the diffusion collision integral

    vector<vector<double> > D_bi;
    D_bi.resize(nsp_);
    for (int isp = 0; isp < nsp_; ++isp)
      D_bi[isp].resize(nsp_, 0.0);

    for ( int isp = 0; isp < nsp_; ++isp ) {
	for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	    if ( isp == jsp ) continue;
            double M_ij = 0.0;
            M_ij = 1/M_[isp] + 1/M_[jsp];
            M_ij = 2/M_ij; //in kg/mol
            M_ij *= 1.0e3; //in g/mol

            double d_c_ij = 0.0;
            d_c_ij = (d_c[isp] + d_c[jsp])/2;

            double e_LJ_ij = 0.0;
            e_LJ_ij = sqrt(e_LJ[isp]*e_LJ[jsp]);

            double T_star = 0.0;
            T_star = Q.T[0]/e_LJ_ij; //dimensionless
            double omega = 0.0;
            omega = Coef_A/(pow(T_star,Coef_B));
            omega += Coef_C/exp(Coef_D*T_star);
            omega += Coef_E/exp(Coef_F*T_star);
            omega += Coef_G/exp(Coef_H*T_star);
            //Diffusion collision integral, dimensionless
            
            D_bi[isp][jsp] = 1.0/(Q.p/1.0e5)/pow(M_ij,1.0/2)/pow(d_c_ij,2.0)/omega;
            //note: pressure unit in bar
            D_bi[isp][jsp] *= 0.00266*pow(Q.T[0],3.0/2);
            //in cm2/s           
            D_bi[isp][jsp] *= 1.0e-4;
            //in m2/s
	    }
	}

    fill_in_x(Q.rho, Q.massf);
    fill_in_DAV_im(D_bi);
    //calculate mixture-averaged diffusion coefficient for species i


    // Set diffusive fluxes...
    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
	jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
	jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }

    if ( FICKS_WITH_SUTTON_AND_GNOFFO_CORRECTION ) {
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
    }

    // Calculate the thermal diffusivity
    //Gas_model *gmodel = get_gas_model_ptr();
    //double k_total = 0.0;
    //for ( size_t itm=0; itm<Q.k.size(); ++itm ) {
    //	k_total += Q.k[itm];
    //}
    //int status;
    //double alpha = k_total / Q.rho / gmodel->Cp(Q, status);

    // Calculate the Lewis number
    //vector<double> Le_DRM19;
    //Le_DRM19.resize(nsp_, 1.0);
    //for ( int isp = 0; isp < nsp_; ++isp ) {
    //    Le_DRM19[isp] = alpha / DAV_im_[isp];
    //}


}


FicksFirstLaw_WD1::
FicksFirstLaw_WD1(const string name, int nsp)
    : DiffusionModel(name, nsp)
{}

FicksFirstLaw_WD1::
FicksFirstLaw_WD1(const FicksFirstLaw_WD1 &f)
    : DiffusionModel(f.name_, f.nsp_)
{}

FicksFirstLaw_WD1::
~FicksFirstLaw_WD1() {}

string
FicksFirstLaw_WD1::
str() const
{
    ostringstream ost;
    cout << "FicksFirstLaw_WD1::str() - not implemented.\n";
    ost << "FicksFirstLaw_WD1";
    return ost.str();
}

void
FicksFirstLaw_WD1::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    //first calculate binary diffusion coefficient D_AB
    //according to the Chapman-Enskog relation
    //refer to "The Properties of Gases and Liquids", Reid et al.

    vector<double> d_c;
    d_c.resize(nsp_, 0.0);
    d_c[0] = 3.621; // N2
    d_c[1] = 3.746; // CH4
    d_c[2] = 3.458; // O2
    d_c[3] = 2.605; // H2O
    d_c[4] = 3.763; // CO2
    //Lennard-Jones collision diameter "sigma" in Angstroms
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    vector<double> e_LJ;
    e_LJ.resize(nsp_, 0.0);
    e_LJ[0] = 97.530; // N2
    e_LJ[1] = 141.400; // CH4
    e_LJ[2] = 107.400; // O2
    e_LJ[3] = 572.400; // H2O
    e_LJ[4] = 244.000; // CO2
    //Lennard-Jones potential well depth "epsilon/kb" in Kelvins
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    double Coef_A = 1.06036;
    double Coef_B = 0.15610;
    double Coef_C = 0.19300;
    double Coef_D = 0.47635;
    double Coef_E = 1.03587;
    double Coef_F = 1.52996;
    double Coef_G = 1.76474;
    double Coef_H = 3.89411;
    //Coefficients needed for calculating the diffusion collision integral

    vector<vector<double> > D_bi;
    D_bi.resize(nsp_);
    for (int isp = 0; isp < nsp_; ++isp)
      D_bi[isp].resize(nsp_, 0.0);

    for ( int isp = 0; isp < nsp_; ++isp ) {
	for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	    if ( isp == jsp ) continue;
            double M_ij = 0.0;
            M_ij = 1/M_[isp] + 1/M_[jsp];
            M_ij = 2/M_ij; //in kg/mol
            M_ij *= 1.0e3; //in g/mol

            double d_c_ij = 0.0;
            d_c_ij = (d_c[isp] + d_c[jsp])/2;

            double e_LJ_ij = 0.0;
            e_LJ_ij = sqrt(e_LJ[isp]*e_LJ[jsp]);

            double T_star = 0.0;
            T_star = Q.T[0]/e_LJ_ij; //dimensionless
            double omega = 0.0;
            omega = Coef_A/(pow(T_star,Coef_B));
            omega += Coef_C/exp(Coef_D*T_star);
            omega += Coef_E/exp(Coef_F*T_star);
            omega += Coef_G/exp(Coef_H*T_star);
            //Diffusion collision integral, dimensionless
            
            D_bi[isp][jsp] = 1.0/(Q.p/1.0e5)/pow(M_ij,1.0/2)/pow(d_c_ij,2.0)/omega;
            //note: pressure unit in bar
            D_bi[isp][jsp] *= 0.00266*pow(Q.T[0],3.0/2);
            //in cm2/s           
            D_bi[isp][jsp] *= 1.0e-4;
            //in m2/s
	    }
	}

    fill_in_x(Q.rho, Q.massf);
    fill_in_DAV_im(D_bi);
    //calculate mixture-averaged diffusion coefficient for species i


    // Set diffusive fluxes...
    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
	jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
	jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }

    if ( FICKS_WITH_SUTTON_AND_GNOFFO_CORRECTION ) {
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
    }

    // Calculate the thermal diffusivity
    //Gas_model *gmodel = get_gas_model_ptr();
    //double k_total = 0.0;
    //for ( size_t itm=0; itm<Q.k.size(); ++itm ) {
    //	k_total += Q.k[itm];
    //}
    //int status;
    //double alpha = k_total / Q.rho / gmodel->Cp(Q, status);

    // Calculate the Lewis number
    //vector<double> Le_WD1;
    //Le_WD1.resize(nsp_, 1.0);
    //for ( int isp = 0; isp < nsp_; ++isp ) {
    //    Le_WD1[isp] = alpha / DAV_im_[isp];
    //}


}


FicksFirstLaw_WD2::
FicksFirstLaw_WD2(const string name, int nsp)
    : DiffusionModel(name, nsp)
{}

FicksFirstLaw_WD2::
FicksFirstLaw_WD2(const FicksFirstLaw_WD2 &f)
    : DiffusionModel(f.name_, f.nsp_)
{}

FicksFirstLaw_WD2::
~FicksFirstLaw_WD2() {}

string
FicksFirstLaw_WD2::
str() const
{
    ostringstream ost;
    cout << "FicksFirstLaw_WD2::str() - not implemented.\n";
    ost << "FicksFirstLaw_WD2";
    return ost.str();
}

void
FicksFirstLaw_WD2::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    //first calculate binary diffusion coefficient D_AB
    //according to the Chapman-Enskog relation
    //refer to "The Properties of Gases and Liquids", Reid et al.

    vector<double> d_c;
    d_c.resize(nsp_, 0.0);
    d_c[0] = 3.621; // N2
    d_c[1] = 3.746; // CH4
    d_c[2] = 3.458; // O2
    d_c[3] = 3.650; // CO
    d_c[4] = 2.605; // H2O
    d_c[5] = 3.763; // CO2
    //Lennard-Jones collision diameter "sigma" in Angstroms
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    vector<double> e_LJ;
    e_LJ.resize(nsp_, 0.0);
    e_LJ[0] = 97.530; // N2
    e_LJ[1] = 141.400; // CH4
    e_LJ[2] = 107.400; // O2
    e_LJ[3] = 98.100; // CO
    e_LJ[4] = 572.400; // H2O
    e_LJ[5] = 244.000; // CO2
    //Lennard-Jones potential well depth "epsilon/kb" in Kelvins
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    double Coef_A = 1.06036;
    double Coef_B = 0.15610;
    double Coef_C = 0.19300;
    double Coef_D = 0.47635;
    double Coef_E = 1.03587;
    double Coef_F = 1.52996;
    double Coef_G = 1.76474;
    double Coef_H = 3.89411;
    //Coefficients needed for calculating the diffusion collision integral

    vector<vector<double> > D_bi;
    D_bi.resize(nsp_);
    for (int isp = 0; isp < nsp_; ++isp)
      D_bi[isp].resize(nsp_, 0.0);

    for ( int isp = 0; isp < nsp_; ++isp ) {
	for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	    if ( isp == jsp ) continue;
            double M_ij = 0.0;
            M_ij = 1/M_[isp] + 1/M_[jsp];
            M_ij = 2/M_ij; //in kg/mol
            M_ij *= 1.0e3; //in g/mol

            double d_c_ij = 0.0;
            d_c_ij = (d_c[isp] + d_c[jsp])/2;

            double e_LJ_ij = 0.0;
            e_LJ_ij = sqrt(e_LJ[isp]*e_LJ[jsp]);

            double T_star = 0.0;
            T_star = Q.T[0]/e_LJ_ij; //dimensionless
            double omega = 0.0;
            omega = Coef_A/(pow(T_star,Coef_B));
            omega += Coef_C/exp(Coef_D*T_star);
            omega += Coef_E/exp(Coef_F*T_star);
            omega += Coef_G/exp(Coef_H*T_star);
            //Diffusion collision integral, dimensionless
            
            D_bi[isp][jsp] = 1.0/(Q.p/1.0e5)/pow(M_ij,1.0/2)/pow(d_c_ij,2.0)/omega;
            //note: pressure unit in bar
            D_bi[isp][jsp] *= 0.00266*pow(Q.T[0],3.0/2);
            //in cm2/s           
            D_bi[isp][jsp] *= 1.0e-4;
            //in m2/s
	    }
	}

    fill_in_x(Q.rho, Q.massf);
    fill_in_DAV_im(D_bi);
    //calculate mixture-averaged diffusion coefficient for species i


    // Set diffusive fluxes...
    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
	jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
	jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }

    if ( FICKS_WITH_SUTTON_AND_GNOFFO_CORRECTION ) {
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
    }

    // Calculate the thermal diffusivity
    //Gas_model *gmodel = get_gas_model_ptr();
    //double k_total = 0.0;
    //for ( size_t itm=0; itm<Q.k.size(); ++itm ) {
    //	k_total += Q.k[itm];
    //}
    //int status;
    //double alpha = k_total / Q.rho / gmodel->Cp(Q, status);

    // Calculate the Lewis number
    //vector<double> Le_WD2;
    //Le_WD2.resize(nsp_, 1.0);
    //for ( int isp = 0; isp < nsp_; ++isp ) {
    //    Le_WD2[isp] = alpha / DAV_im_[isp];
    //}


}

FicksFirstLaw_JL4::
FicksFirstLaw_JL4(const string name, int nsp)
    : DiffusionModel(name, nsp)
{}

FicksFirstLaw_JL4::
FicksFirstLaw_JL4(const FicksFirstLaw_JL4 &f)
    : DiffusionModel(f.name_, f.nsp_)
{}

FicksFirstLaw_JL4::
~FicksFirstLaw_JL4() {}

string
FicksFirstLaw_JL4::
str() const
{
    ostringstream ost;
    cout << "FicksFirstLaw_JL4::str() - not implemented.\n";
    ost << "FicksFirstLaw_JL4";
    return ost.str();
}

void
FicksFirstLaw_JL4::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    //first calculate binary diffusion coefficient D_AB
    //according to the Chapman-Enskog relation
    //refer to "The Properties of Gases and Liquids", Reid et al.

    vector<double> d_c;
    d_c.resize(nsp_, 0.0);
    d_c[0] = 3.621; // N2
    d_c[1] = 3.746; // CH4
    d_c[2] = 3.458; // O2
    d_c[3] = 3.650; // CO
    d_c[4] = 2.920; // H2
    d_c[5] = 2.605; // H2O
    d_c[6] = 3.763; // CO2
    //Lennard-Jones collision diameter "sigma" in Angstroms
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    vector<double> e_LJ;
    e_LJ.resize(nsp_, 0.0);
    e_LJ[0] = 97.530; // N2
    e_LJ[1] = 141.400; // CH4
    e_LJ[2] = 107.400; // O2
    e_LJ[3] = 98.100; // CO
    e_LJ[4] = 38.000; // H2
    e_LJ[5] = 572.400; // H2O
    e_LJ[6] = 244.000; // CO2
    //Lennard-Jones potential well depth "epsilon/kb" in Kelvins
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    double Coef_A = 1.06036;
    double Coef_B = 0.15610;
    double Coef_C = 0.19300;
    double Coef_D = 0.47635;
    double Coef_E = 1.03587;
    double Coef_F = 1.52996;
    double Coef_G = 1.76474;
    double Coef_H = 3.89411;
    //Coefficients needed for calculating the diffusion collision integral

    vector<vector<double> > D_bi;
    D_bi.resize(nsp_);
    for (int isp = 0; isp < nsp_; ++isp)
      D_bi[isp].resize(nsp_, 0.0);

    for ( int isp = 0; isp < nsp_; ++isp ) {
	for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	    if ( isp == jsp ) continue;
            double M_ij = 0.0;
            M_ij = 1/M_[isp] + 1/M_[jsp];
            M_ij = 2/M_ij; //in kg/mol
            M_ij *= 1.0e3; //in g/mol

            double d_c_ij = 0.0;
            d_c_ij = (d_c[isp] + d_c[jsp])/2;

            double e_LJ_ij = 0.0;
            e_LJ_ij = sqrt(e_LJ[isp]*e_LJ[jsp]);

            double T_star = 0.0;
            T_star = Q.T[0]/e_LJ_ij; //dimensionless
            double omega = 0.0;
            omega = Coef_A/(pow(T_star,Coef_B));
            omega += Coef_C/exp(Coef_D*T_star);
            omega += Coef_E/exp(Coef_F*T_star);
            omega += Coef_G/exp(Coef_H*T_star);
            //Diffusion collision integral, dimensionless
            
            D_bi[isp][jsp] = 1.0/(Q.p/1.0e5)/pow(M_ij,1.0/2)/pow(d_c_ij,2.0)/omega;
            //note: pressure unit in bar
            D_bi[isp][jsp] *= 0.00266*pow(Q.T[0],3.0/2);
            //in cm2/s           
            D_bi[isp][jsp] *= 1.0e-4;
            //in m2/s
	    }
	}

    fill_in_x(Q.rho, Q.massf);
    fill_in_DAV_im(D_bi);
    //calculate mixture-averaged diffusion coefficient for species i


    // Set diffusive fluxes...
    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
	jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
	jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }

    if ( FICKS_WITH_SUTTON_AND_GNOFFO_CORRECTION ) {
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
    }

    // Calculate the thermal diffusivity
    //Gas_model *gmodel = get_gas_model_ptr();
    //double k_total = 0.0;
    //for ( size_t itm=0; itm<Q.k.size(); ++itm ) {
    //	k_total += Q.k[itm];
    //}
    //int status;
    //double alpha = k_total / Q.rho / gmodel->Cp(Q, status);

    // Calculate the Lewis number
    //vector<double> Le_JL4;
    //Le_JL4.resize(nsp_, 1.0);
    //for ( int isp = 0; isp < nsp_; ++isp ) {
    //    Le_JL4[isp] = alpha / DAV_im_[isp];
    //}


}


FicksFirstLaw_YSSS5::
FicksFirstLaw_YSSS5(const string name, int nsp)
    : DiffusionModel(name, nsp)
{}

FicksFirstLaw_YSSS5::
FicksFirstLaw_YSSS5(const FicksFirstLaw_YSSS5 &f)
    : DiffusionModel(f.name_, f.nsp_)
{}

FicksFirstLaw_YSSS5::
~FicksFirstLaw_YSSS5() {}

string
FicksFirstLaw_YSSS5::
str() const
{
    ostringstream ost;
    cout << "FicksFirstLaw_YSSS5::str() - not implemented.\n";
    ost << "FicksFirstLaw_YSSS5";
    return ost.str();
}

void
FicksFirstLaw_YSSS5::
calculate_diffusion_fluxes(const Gas_data &Q,
			   double D_t,
			   const vector<double> &dfdx, 
			   const vector<double> &dfdy,
			   const vector<double> &dfdz,
			   vector<double> &jx, 
			   vector<double> &jy,
			   vector<double> &jz)
{
    //first calculate binary diffusion coefficient D_AB
    //according to the Chapman-Enskog relation
    //refer to "The Properties of Gases and Liquids", Reid et al.

    vector<double> d_c;
    d_c.resize(nsp_, 0.0);
    d_c[0] = 3.621; // N2
    d_c[1] = 3.626; // CH3OH
    d_c[2] = 3.458; // O2
    d_c[3] = 3.763; // CO2
    d_c[4] = 2.605; // H2O
    d_c[5] = 3.590; // CH2O
    d_c[6] = 3.650; // CO
    d_c[7] = 2.920; // H2
    d_c[8] = 2.050; // H
    //Lennard-Jones collision diameter "sigma" in Angstroms
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    vector<double> e_LJ;
    e_LJ.resize(nsp_, 0.0);
    e_LJ[0] = 97.530; // N2
    e_LJ[1] = 481.800; // CH3OH
    e_LJ[2] = 107.400; // O2
    e_LJ[3] = 244.000; // CO2
    e_LJ[4] = 572.400; // H2O
    e_LJ[5] = 498.000; // CH2O
    e_LJ[6] = 98.100; // CO
    e_LJ[7] = 38.000; // H2
    e_LJ[8] = 145.000; // H
    //Lennard-Jones potential well depth "epsilon/kb" in Kelvins
    //refer to GRI-Mech 3.0 transport file (Chemkin format)
    //http://combustion.berkeley.edu/gri-mech/version30/text30.html

    double Coef_A = 1.06036;
    double Coef_B = 0.15610;
    double Coef_C = 0.19300;
    double Coef_D = 0.47635;
    double Coef_E = 1.03587;
    double Coef_F = 1.52996;
    double Coef_G = 1.76474;
    double Coef_H = 3.89411;
    //Coefficients needed for calculating the diffusion collision integral

    vector<vector<double> > D_bi;
    D_bi.resize(nsp_);
    for (int isp = 0; isp < nsp_; ++isp)
      D_bi[isp].resize(nsp_, 0.0);

    for ( int isp = 0; isp < nsp_; ++isp ) {
	for ( int jsp = 0; jsp < nsp_; ++jsp ) {
	    if ( isp == jsp ) continue;
            double M_ij = 0.0;
            M_ij = 1/M_[isp] + 1/M_[jsp];
            M_ij = 2/M_ij; //in kg/mol
            M_ij *= 1.0e3; //in g/mol

            double d_c_ij = 0.0;
            d_c_ij = (d_c[isp] + d_c[jsp])/2;

            double e_LJ_ij = 0.0;
            e_LJ_ij = sqrt(e_LJ[isp]*e_LJ[jsp]);

            double T_star = 0.0;
            T_star = Q.T[0]/e_LJ_ij; //dimensionless
            double omega = 0.0;
            omega = Coef_A/(pow(T_star,Coef_B));
            omega += Coef_C/exp(Coef_D*T_star);
            omega += Coef_E/exp(Coef_F*T_star);
            omega += Coef_G/exp(Coef_H*T_star);
            //Diffusion collision integral, dimensionless
            
            D_bi[isp][jsp] = 1.0/(Q.p/1.0e5)/pow(M_ij,1.0/2)/pow(d_c_ij,2.0)/omega;
            //note: pressure unit in bar
            D_bi[isp][jsp] *= 0.00266*pow(Q.T[0],3.0/2);
            //in cm2/s           
            D_bi[isp][jsp] *= 1.0e-4;
            //in m2/s
	    }
	}

    fill_in_x(Q.rho, Q.massf);
    fill_in_DAV_im(D_bi);
    //calculate mixture-averaged diffusion coefficient for species i


    // Set diffusive fluxes...
    for ( int isp = 0; isp < nsp_; ++isp ) {
	jx[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdx[isp];
	jy[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdy[isp];
	jz[isp] = -Q.rho * (DAV_im_[isp] + D_t) * dfdz[isp];
    }

    if ( FICKS_WITH_SUTTON_AND_GNOFFO_CORRECTION ) {
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
    }

    // Calculate the thermal diffusivity
    //Gas_model *gmodel = get_gas_model_ptr();
    //double k_total = 0.0;
    //for ( size_t itm=0; itm<Q.k.size(); ++itm ) {
    //	k_total += Q.k[itm];
    //}
    //int status;
    //double alpha = k_total / Q.rho / gmodel->Cp(Q, status);

    // Calculate the Lewis number
    //vector<double> Le_YSSS5;
    //Le_YSSS5.resize(nsp_, 1.0);
    //for ( int isp = 0; isp < nsp_; ++isp ) {
    //    Le_YSSS5[isp] = alpha / DAV_im_[isp];
    //}


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
	global_data *gd = get_global_data_ptr();
	dmodel = new ConstantLewisNumber("Constant Lewis number diffusion", nsp, gd->diffusion_lewis);
    }
    else if ( diffusion_model == "ConstantSchmidtNumber" ) {
	global_data *gd = get_global_data_ptr();
        dmodel = new ConstantSchmidtNumber("Constant Schmidt number diffusion", nsp, gd->diffusion_schmidt);
    }
    else if ( diffusion_model == "ConstantLewisNumber_DRM19" ) {
	global_data *gd = get_global_data_ptr();
	dmodel = new ConstantLewisNumber_DRM19("Constant Lewis number diffusion-DRM19 scheme", nsp );
    }
    else if ( diffusion_model == "FicksFirstLaw_DRM19" ) {
	dmodel = new FicksFirstLaw_DRM19("Fick's first law (of diffusion) for DRM19 scheme", nsp );
    }
    else if ( diffusion_model == "FicksFirstLaw_WD1" ) {
	dmodel = new FicksFirstLaw_WD1("Fick's first law (of diffusion) for WD1 scheme", nsp );
    }
    else if ( diffusion_model == "FicksFirstLaw_WD2" ) {
	dmodel = new FicksFirstLaw_WD2("Fick's first law (of diffusion) for WD2 scheme", nsp );
    }
    else if ( diffusion_model == "FicksFirstLaw_JL4" ) {
	dmodel = new FicksFirstLaw_JL4("Fick's first law (of diffusion) for JL4 scheme", nsp );
    }
    else if ( diffusion_model == "FicksFirstLaw_YSSS5" ) {
	dmodel = new FicksFirstLaw_YSSS5("Fick's first law (of diffusion) for YSSS5 scheme", nsp );
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

