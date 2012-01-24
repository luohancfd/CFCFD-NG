// Author: Daniel F Potter
// Version: 12-Oct-2009
//          Initial coding.
// Version: 13-Apr-2010
//          Changed name from CI_functor to GuptaYos_CI_functor.

#include <iostream>
#include <sstream>
#include <cmath>

#include "../../util/source/lua_service.hh"

#include "CI-functor.hh"
#include "physical_constants.hh"

using namespace std;

/****************************** Gupta-Yos curve fits **************************/

GuptaYos_CI_functor::
GuptaYos_CI_functor( string section, double T_low, double T_high, lua_State *L )
 : T_low_( T_low ), T_high_( T_high )
{
    if ( T_high_ <= T_low_ ) {
    	ostringstream ost;
	ost << "GuptaYos_CI_functor::GuptaYos_CI_functor()\n";
	ost << "Error initialising GuptaYos_CI_functor object.\n";
	ost << "T_high should be greater than T_low\n";
	ost << "T_high= " << T_high_ << " T_low= " << T_low_ << endl;
	ost << "Bailing out!\n";
	input_error(ost);
    }

    lua_getfield(L, -1, section.c_str() );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "GuptaYos_CI_functor::GuptaYos_CI_functor()\n";
	ost << "A '" << section << "' table is expected, but it's not found.\n";
	input_error(ost);
    }

    int n = lua_objlen(L, -1);
    if ( n != 4 ) {
    	ostringstream ost;
	ost << "GuptaYos_CI_functor::GuptaYos_CI_functor()\n";
	ost << "4 coefficients are required to specify the\n";
	ost << "collision integral curves in the 'coeffs' table.\n";
	ost << "Bailing out!\n";
	input_error(ost);
    }
    
    a_.resize(4);
    
    for ( size_t i = 0; i < a_.size(); ++i ) {
	lua_rawgeti(L, -1, i+1);
	a_[i] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }

    lua_pop(L, 1); // pops 'section'
}

GuptaYos_CI_functor::
GuptaYos_CI_functor(const GuptaYos_CI_functor &c)
    : T_low_(c.T_low_), T_high_(c.T_high_), a_(c.a_)  {}
    
GuptaYos_CI_functor::
~GuptaYos_CI_functor() {}

GuptaYos_CI_functor*
GuptaYos_CI_functor::
clone() const
{
    return new GuptaYos_CI_functor(*this);
}


double
GuptaYos_CI_functor::
operator()(double T)
{
    // 0. Impose temperature limits
    if ( T < T_low_ ) T = T_low_;
    if ( T > T_high_ ) T = T_high_;
    
    // 1. Evaluate curve at T
    return exp( a_[3] ) * pow( T, a_[0]*pow( log(T),2) + a_[1]*log(T) + a_[2] );
}

/************************* Bruno's curve fits **************************/

/* Heavy particle interactions (at least one neutral species) */

Bruno_HPI_CI_functor::
Bruno_HPI_CI_functor( int l, int s, int Z_l, int Z_s, lua_State * L )
{
    eps0_ = get_positive_number(L,-1,"eps0") * 0.001 * PC_e_SI;	// convert meV -> J
    double beta = get_positive_number(L,-1,"beta");
    double re = get_positive_number(L,-1,"re");
    
    double x0 = 0.0;
    if ( Z_l==0 && Z_s==0 ) {
    	// Neutral-neutral interaction
    	x0 = 0.8002 * pow( beta, 0.049256 );
    }
    else if ( ( Z_l==1 && Z_s==0 ) || ( Z_l==0 && Z_s==1 ) ) {
    	// Neutral-ion interaction
    	x0 = 0.7564 * pow( beta, 0.064605 );
    }
    else {
    	ostringstream oss;
    	oss << "Bruno_HPI_CI_functor::Bruno_HPI_CI_functor()" << endl
    	    << "Expected a neutral-neutral or neutral-ion interaction" << endl
    	    << "Z_l = " << Z_l << ", Z_s = " << Z_s << endl;
    	input_error(oss);
    }
    sigma2_ = pow( x0*re, 2 );		// Collision diameter in Ang**2

    a_.resize(7);
    vector<double> c0(7), c1(7), c2(7);
    if ( l==1 && s==1 ) {
    	c0[0] =  7.884756e-1; c1[0] = -2.438494e-2; c2[0] = 0.0;
    	c0[1] = -2.952759e-1; c1[1] = -1.744149e-3; c2[1] = 0.0;
    	c0[2] =  5.020892e-1; c1[2] =  4.316985e-2; c2[2] = 0.0;
    	c0[3] = -9.042460e-1; c1[3] = -4.017103e-2; c2[3] = 0.0;
    	c0[4] = -3.373058e+0; c1[4] =  2.458535e-1; c2[4] = -4.850047e-3;
    	c0[5] =  4.161981e+0; c1[5] =  2.202737e-1; c2[5] = -1.718010e-2;
    	c0[6] =  2.462523e+0; c1[6] =  3.231308e-1; c2[6] = -2.281072e-2;
    }
    else if ( l==2 && s==2 ) {
    	c0[0] =  7.898524e-1; c1[0] = -2.114115e-2; c2[0] = 0.0;
    	c0[1] = -2.998325e-1; c1[1] = -1.243977e-3; c2[1] = 0.0;
    	c0[2] =  7.077103e-1; c1[2] =  3.583907e-2; c2[2] = 0.0;
    	c0[3] = -8.946857e-1; c1[3] = -2.473947e-2; c2[3] = 0.0;
    	c0[4] = -2.958969e+0; c1[4] =  2.303358e-1; c2[4] = -5.226562e-3;
    	c0[5] =  4.348412e+0; c1[5] =  1.920321e-1; c2[5] = -1.496557e-2;
    	c0[6] =  2.205440e+0; c1[6] =  2.567027e-1; c2[6] = -1.861359e-2;
    }
    else {
    	ostringstream oss;
    	oss << "Bruno_NHPI_CI_functor::Bruno_NHPI_CI_functor()" << endl
    	    << "l-s combination: " << l << "-" << s << " not available." << endl;
    	input_error(oss);
    }
    for ( int i=0; i<7; ++i )
    	a_[i] = c0[i] + c1[i]*beta + c2[i]*beta;
}

Bruno_HPI_CI_functor::
Bruno_HPI_CI_functor(const Bruno_HPI_CI_functor &c)
    : eps0_( c.eps0_ ), a_( c.a_ ) {}
    
Bruno_HPI_CI_functor::
~Bruno_HPI_CI_functor() {}

Bruno_HPI_CI_functor*
Bruno_HPI_CI_functor::
clone() const
{
    return new Bruno_HPI_CI_functor(*this);
}

double
Bruno_HPI_CI_functor::
operator()(double T)
{
    // Eq. 5.4 with sigma**2 at the front to get the dimensional cross-section
    double x = log( PC_k_SI * T / eps0_ );
    double tmpA = a_[0] + a_[1]*x;
    double tmpB = exp( ( x - a_[2] ) / a_[3] );
    double tmpC = exp( ( a_[2] - x ) / a_[3] );
    double tmpD = exp( ( x - a_[5] ) / a_[6] );
    double tmpE = exp( ( a_[5] - x ) / a_[6] );
    
    return M_PI * sigma2_ * exp( tmpA * tmpB / ( tmpB + tmpC ) + a_[4] * tmpD / ( tmpD + tmpE ) );
}

/* Charged-charged particle interaction */

Bruno_CPI_CI_functor::
Bruno_CPI_CI_functor( int l, int m, lua_State * L )
{
    // NOTE: eps0 is not given in the Bruno paper...
    eps0_ = get_positive_number(L,-1,"eps0") * 0.001 * PC_e_SI;	// convert meV -> J
    
    double Nl = 1.0 / ( 1.0 - ( 1.0 + pow( -1.0, l ) ) / ( 2.0 * ( double(l) + 1.0 ) ) );
    constA_ = double(l) * Nl / ( double(m) * ( double(m) + 1.0 ) );
    
    double Am = 0.0;
    if ( m > 1 ) {
    	for ( int i=2; i<=m; ++i ) Am += 1.0 / ( double(i) - 1.0 );
    }
    
    double Cl = 0.0;
    if ( l%2 ) {
    	// l is odd
    	for ( int i=1; i<=l; i+=2 ) Cl+=1.0/double(i);
    	Cl -= 1.0 / ( 2.0*double(l) );
    }
    else {
    	// l is even
    	for ( int i=2; i<=l; i+=2 ) Cl+=1.0/(double(i)-1.0);
    }
    
    constB_ = 4.0 / ( exp( 0.57721 ) * exp( 0.57721 ) ) * exp( Am - Cl ); 
}

Bruno_CPI_CI_functor::
Bruno_CPI_CI_functor(const Bruno_CPI_CI_functor &c)
    : constA_(c.constA_), constB_(c.constB_) {}
    
Bruno_CPI_CI_functor::
~Bruno_CPI_CI_functor() {}

Bruno_CPI_CI_functor*
Bruno_CPI_CI_functor::
clone() const
{
    return new Bruno_CPI_CI_functor(*this);
}

double
Bruno_CPI_CI_functor::
operator()(double T, double N_e)
{
    // Eq. 5.25
    // Firstly calculate the Debye length (using wrights formula for the moment)
    double lambda_D_CGS = sqrt( PC_k_CGS * T / ( 4.0 * M_PI * N_e * 1.0e-6 * PC_e_CGS * PC_e_CGS ) );
    double sigma = 2.0 * lambda_D_CGS * 1.0e-8;		// convert cm -> Ang
    
    double T_star = T * PC_k_SI / eps0_;
    
    return sigma * sigma * constA_ * log( constB_ * T_star + 1.0 ) / T_star / T_star;
}

/***************************** Wright curve fits ****************************/

Stallcop_SCP_CI_functor::
Stallcop_SCP_CI_functor( int l, int m, int Z_l, int Z_m )
{
    if ( l==1 && m==1 ) {
    	// Diffusion CI
    	if ( Z_l==Z_m ) {
    	    // Repulsive potential
    	    C_ = 0.138; c_ = 0.0106; D_ = 0.765;
    	}
    	else {
    	    // attractive potential
    	    C_ = -0.476; c_ = 0.0313; D_ = 0.784;
    	}
    }
    else if ( l==2 && m==2 ) {
    	// Viscosity CI
    	if ( Z_l==Z_m ) {
    	    // Repulsive potential
    	    C_ = 0.157; c_ = 0.0274; D_ = 1.235;
    	}
    	else {
    	    // attractive potential
    	    C_ = -0.146; c_ = 0.0377; D_ = 1.262;
    	}
    }
}

Stallcop_SCP_CI_functor::
Stallcop_SCP_CI_functor(const Stallcop_SCP_CI_functor &c)
    : C_(c.C_), c_(c.c_), D_(c.D_) {}
    
Stallcop_SCP_CI_functor::
~Stallcop_SCP_CI_functor() {}

Stallcop_SCP_CI_functor*
Stallcop_SCP_CI_functor::
clone() const
{
    return new Stallcop_SCP_CI_functor(*this);
}

double
Stallcop_SCP_CI_functor::
operator()(double T, double N_e)
{
    // Eq. 6
    // NOTE: - think we have to use CGS units for this one
    // Firstly calculate the Debye length
    double lambda_D = sqrt( PC_k_CGS * T / ( 4.0 * M_PI * N_e * 1.0e-6 * PC_e_CGS * PC_e_CGS ) );	// cm
    // double lambda_D = 6.905*sqrt(T/(N_e*1.0e-6));
    // cout << "lambda_D = " << lambda_D << endl;
    
    // Now the reduced temperature
    double T_star = lambda_D / ( PC_e_CGS*PC_e_CGS/PC_k_CGS/T );					// ND
    // double T_star = 4132.5*pow(T,1.5)/sqrt(N_e*1.0e-6);
    // cout << "T_star = " << T_star << endl;
    
    return M_PI*5.0e15*(lambda_D/T_star)*(lambda_D/T_star)*log(D_*T_star*(1.0-C_*exp(-c_*T_star))+1.0);	// Ang**2
}

/************************* Neufeld's curve fits **************************/

Neufeld_CI_functor::
Neufeld_CI_functor( int l, int s, Chemical_species * I, Chemical_species * J, lua_State * L )
{
    // Leonard Jones potential parameters
    eps0_K_ = sqrt( I->get_eps0() * J->get_eps0() ) / PC_k_SI;			// Potential well in kelvin
    sigma_ = 0.5 * ( I->get_sigma() + J->get_sigma() ) * 1.0e+10;		// Rigid sphere collision diameter in angstroms

    if ( eps0_K_==0.0 || sigma_==0.0 ) {
    	ostringstream oss;
    	oss << "Neufeld_CI_functor::Neufeld_CI_functor" << endl
    	    << "eps0_K_ = " << eps0_K_ << ", sigma_ = " << sigma_ << endl
    	    << "The LJ parameters for the colliding species ("
    	    << I->get_name() << ", " << J->get_name() << ") are not in the species library" << endl;
    	input_error(oss);
    }

    c_.resize(12);
    if ( l==1 && s==1 ) {
    	c_[0] = 1.06036;   c_[1] = 0.15610; c_[2]  = 0.19300;  c_[3]  = 0.47635;
    	c_[4] = 1.03587;   c_[5] = 1.52996; c_[6]  = 1.76474;  c_[7]  = 3.89411;
    	c_[8] = 0.00000;   c_[9] = 0.00000; c_[10] = 0.00000;  c_[11] = 0.00000;
    }
    else if ( l==2 && s==2 ) {
    	c_[0] =  1.16145;  c_[1] = 0.14874; c_[2]  = 0.52487;  c_[3]  = 0.77320;
    	c_[4] =  2.16178;  c_[5] = 2.43787; c_[6]  = 0.00000;  c_[7]  = 0.00000;
    	c_[8] = -6.435e-4; c_[9] = 18.0323; c_[10] = -0.76830; c_[11] = 7.27371;
    }
    else {
    	ostringstream oss;
    	oss << "Neufeld_CI_functor::Neufeld_CI_functor()" << endl
    	    << "l-s combination: " << l << "-" << s << " not available." << endl;
    	input_error(oss);
    }
}

Neufeld_CI_functor::
Neufeld_CI_functor(const Neufeld_CI_functor &n)
    : eps0_K_( n.eps0_K_ ), sigma_( n.sigma_ ), c_( n.c_ ) {}
    
Neufeld_CI_functor::
~Neufeld_CI_functor() {}

Neufeld_CI_functor*
Neufeld_CI_functor::
clone() const
{
    return new Neufeld_CI_functor(*this);
}

double
Neufeld_CI_functor::
operator()(double T)
{
    // See Table I of Neufeld et al (1972) JCP 57(3)
    
    // 0. Reduced temperature
    double T_star = T / ( eps0_K_ );
    
    // 1. Impose reduced temperature range
    // if ( T_star < 0.3 ) T_star = 0.3;
    // else if ( T_star > 100.0 ) T_star = 100.0;
    
    // 2. Calculate the reduced collision integral
    double omega_star = c_[0] / pow( T_star, c_[1] ) + c_[2] / exp( c_[3] * T_star ) + \
    			c_[4] / exp( c_[5] * T_star ) + c_[6] / exp( c_[7] * T_star ) + \
    			c_[8] * pow( T_star, c_[1] ) * sin( c_[9] * pow( T_star, c_[10] ) - c_[11] );
    
    // 3. Return pi x sigma**2 x omega_star in Angstroms**2
    return M_PI * sigma_ * sigma_ * omega_star;
}
