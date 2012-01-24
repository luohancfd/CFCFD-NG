/*  \file Rowan J Gollan
 *  \brief Definitions for the RateCoeffModels (and derivatives)
 *
 *  \author Rowan J Gollan
 *  \version 24-Feb-2006
 **/

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include "../../util/source/config_parser.hh"
#include "../../util/source/useful.h"
#include "../models/gas.hh"
#include "../models/physical_constants.hh"
#include "rr_coeffs.hh"
using namespace std;

RateCoeffModel::RateCoeffModel()
    : model_( string("") ), suffix_( string("") ), k_( 0.0 ) {}

RateCoeffModel::RateCoeffModel( const string model, const string suffix )
    : model_( model ), suffix_( suffix ), k_( 0.0 ) {}

RateCoeffModel::RateCoeffModel( const RateCoeffModel &r )
    : model_( r.model_ ), suffix_( r.suffix_ ), k_( 0.0 ) {}

RateCoeffModel::~RateCoeffModel() {}

string RateCoeffModel::str() const
{
    ostringstream ost;
    ost << "";
    return ost.str();

}

RateCoeffModel* RateCoeffModel::clone()
{
    return new RateCoeffModel( *this );
}

void RateCoeffModel::eval( const Gas_data &Q )
{
    cout << "RateCoeffModel::eval() ---\n";
    cout << "This function is a dummy only! It should not be called ever.\n";
    return;
    
}

GeneralisedArrhenius::
GeneralisedArrhenius()
    : RateCoeffModel(),
      A_(0.0), n_(0.0), E_a_(0.0) {}

GeneralisedArrhenius::
GeneralisedArrhenius( const string suffix, double A, double n, double E_a )
    : RateCoeffModel( "GeneralisedArrhenius", suffix),
      A_( A ), n_( n ), E_a_( E_a ) {}

GeneralisedArrhenius::
GeneralisedArrhenius( ConfigParser &cfg, const string section,
		      const string suffix )
    : RateCoeffModel( "GeneralisedArrhenius", suffix)
{
    string k = "k" + suffix;
    vector<double> val;
    vector<double> notfound;
    notfound.push_back(-1.0);
    if( ! cfg.parse_vector_of_doubles( section, k, val, notfound ) ) {
	cerr << "GeneralisedArrhenius::GeneralisedArrhenius() --- "
	     <<  __FILE__ << (__LINE__ - 2) << endl
	     << "Error reading " << k << " in [" << section << "] of: "
	     << cfg.file_name
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    // Check the value makes sense.
    if( val.size() == 3 ) { // That's a good start
	A_ = val[0];
	n_ = val[1];
	E_a_ = val[2];
    }
    else {
	cerr << "GeneralisedArrhenius::GeneralisedArrhenius() --- "
	     <<  __FILE__ << (__LINE__ - 8) << endl
	     << "Error with number of values listed for GeneralisedArrhenius\n"
	     << "coefficient parameters.\n"
	     << "Three values are exepected: A  n  E_a in expression A T^n exp(-E_a/(kT))\n"
	     << "A size of value vector " << val.size() << " has been found.\n"
	     << "Perhaps the key " << k << " has not been listed in [" << section << "] "
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    return;
}

GeneralisedArrhenius::
GeneralisedArrhenius( const GeneralisedArrhenius &g )
    : RateCoeffModel( g.model_, g.suffix_ ),
      A_( g.A_ ), n_( g.n_ ), E_a_( g.E_a_ ) {}

GeneralisedArrhenius::
~GeneralisedArrhenius() {}

GeneralisedArrhenius*
GeneralisedArrhenius::
clone()
{
    return new GeneralisedArrhenius(*this);
}

GeneralisedArrhenius& 
GeneralisedArrhenius::operator=(const GeneralisedArrhenius &g) 
{
    if (this == &g) // avoid aliasing
	return *this; 
    else {
	A_ = g.A_;
	n_ = g.n_;
	E_a_ = g.E_a_;
    }
    return (*this);
}


string
GeneralisedArrhenius::
str() const
{
    ostringstream ost;
    ost << setprecision(6) << showpoint;
    ost << "k" << suffix_ << "                        = "
	<< A_ << "  "
	<< n_ << "  "
	<< E_a_ <<  endl;
    return ost.str();
}

string
GeneralisedArrhenius::
latex_str(int nc) const
{
    const double base = 1.0e6; // for converting mks -> cgs
    double conv_factor = pow( base, nc-1 );
    double A = A_*conv_factor;

    double E = E_a_/PC_k_SI;

    ostringstream ost;
    ost << showpoint;
    ost << "$k" << suffix_ << "$ & "
	<< setprecision(3) << scientific << A << " & "
	<< setprecision(2) << fixed << n_ << " & "
	<< setprecision(1) << E;
    return ost.str();
}


void
GeneralisedArrhenius::
eval( const Gas_data &Q )
{
    k_ = A_ * pow(Q.T, n_) * exp(- E_a_ / (PC_k_SI * Q.T) );
}



HeavisideZeroActivationEnergy::
HeavisideZeroActivationEnergy( const string suffix, double alpha, double T_i )
    : RateCoeffModel( "HeavisideZeroActivationEnergy", suffix),
      alpha_( alpha ), T_i_( T_i ) {}

HeavisideZeroActivationEnergy::
HeavisideZeroActivationEnergy( ConfigParser &cfg, const string section,
		      const string suffix )
    : RateCoeffModel( "HeavisideZeroActivationEnergy", suffix)
{
    string k = "k" + suffix;
    vector<double> val;
    vector<double> notfound;
    notfound.push_back(-1.0);
    if( ! cfg.parse_vector_of_doubles( section, k, val, notfound ) ) {
	cerr << "HeavisideZeroActivationEnergy::HeavisideZeroActivationEnergy() --- "
	     << "Error reading " << k << " in [" << section << "] of: "
	     << cfg.file_name
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    // Check the value makes sense.
    if( val.size() == 2 ) { // That's a good start
	alpha_ = val[0];
	T_i_ = val[1];
    }
    else {
	cerr << "HeavisideZeroActivationEnergy::HeavisideZeroActivationEnergy() --- "
	     << "Error with number of values listed for HeavisideZeroActivationEnergy\n"
	     << "coefficient parameters.\n"
	     << "Two values are exepected: alpha and T_i in expression alpha * c[i] * H(T - T_i)\n"
	     << "A size of value vector " << val.size() << " has been found.\n"
	     << "Perhaps the key " << k << " has not been listed in [" << section << "] "
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    return;
}

HeavisideZeroActivationEnergy::
HeavisideZeroActivationEnergy( const HeavisideZeroActivationEnergy &g )
    : RateCoeffModel( g.model_, g.suffix_ ),
      alpha_( g.alpha_ ), T_i_( g.T_i_ ) {}

HeavisideZeroActivationEnergy::
~HeavisideZeroActivationEnergy() {}

HeavisideZeroActivationEnergy*
HeavisideZeroActivationEnergy::
clone()
{
    return new HeavisideZeroActivationEnergy(*this);
}

string
HeavisideZeroActivationEnergy::
str() const
{
    ostringstream ost;
    ost << setprecision(6) << showpoint;
    ost << "k" << suffix_ << "                        = "
	<< alpha_ << "  "
	<< T_i_ <<  endl;
    return ost.str();
}

void
HeavisideZeroActivationEnergy::
eval( const Gas_data &Q )
{
    // The if-else statement represents the Heaviside step function
    if( Q.T >= T_i_ ) {
	k_ = alpha_;
    }
    else {
	k_ = 0.0;
    }
}

Nonequilibrium_dissociation::
Nonequilibrium_dissociation(const string suffix,
			    double A, double n, double E_a,
			    int ivib, double alpha, double U)
    : GeneralisedArrhenius(suffix, A, n, E_a),
      ivib_(ivib), alpha_(alpha), U_(U), D_(E_a)
{
    model_ = "Nonequilibrium_dissociation";

}

Nonequilibrium_dissociation::
Nonequilibrium_dissociation(ConfigParser &cfg,
			    const string section, const string suffix)
    : GeneralisedArrhenius(cfg, section, suffix)
{
    D_ = E_a_;
    model_ = "Nonequilibrium_dissociation";
    string alpha_s = "alpha" + suffix;
    if( ! cfg.parse_double(section, alpha_s, alpha_, 0.0) ) {
	cout << "Error reading in value " << alpha_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string ivib_s = "ivib" + suffix;
    

    if( ! cfg.parse_int(section, ivib_s, ivib_, -1) ) {
	cout << "Error reading in value for " << ivib_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( ivib_ < 0 ) {
	cout << "Error reading in value for ivib_f in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "The value should be listed as an integer, greater than or equal to 0.\n"
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string U_s = "U" + suffix;

    if( ! cfg.parse_double(section, U_s, U_, 0.0) ) {
	cout << "Error reading in value " << U_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

}

Nonequilibrium_dissociation::
Nonequilibrium_dissociation(const Nonequilibrium_dissociation &n)
    :  GeneralisedArrhenius(n.suffix_, n.A_, n.n_, n.E_a_),
       ivib_(n.ivib_), alpha_(n.alpha_), U_(n.U_), D_(n.E_a_)
{
    model_ = "Nonequilibrium_dissociation";

}

Nonequilibrium_dissociation::
~Nonequilibrium_dissociation() {}


Nonequilibrium_dissociation*
Nonequilibrium_dissociation::
clone()
{
    return new Nonequilibrium_dissociation(*this);
}

string
Nonequilibrium_dissociation::
str() const
{
    ostringstream ost;
    ost << GeneralisedArrhenius::str()
	<< "alpha" << suffix_ << "                    = " << alpha_ << endl
	<< "ivib" << suffix_ << "                     = " << ivib_ << endl
	<< "U" << suffix_ << "                        = " << U_ << endl
	<< "D" << suffix_ << "                        = " << E_a_ << endl
	<< "A" << suffix_ << "                        = " << E_a_ << endl;
    return ost.str();
}

void
Nonequilibrium_dissociation::
eval(const Gas_data &Q)
{

    // This call set k_
    GeneralisedArrhenius::eval(Q);
    double psi;
    if ( polyatomic_check(ivib_) ) {
        valarray<double> theta_vs;
        get_theta_vs(theta_vs,ivib_);
        psi = polyatomic_nonequilibrium_factor(Q, ivib_, theta_vs, alpha_, D_, U_, D_);
    } else
        psi = nonequilibrium_factor(Q, ivib_, get_theta_v(ivib_), alpha_, D_, U_, D_);
    // cout << "Nonequilibrium_dissociation::eval()\n";
//     cout << "k_eq= " << k_ << endl;
    k_ = psi * k_;
    // cout << "psi= " << psi << " k= " << k_ << endl;

}


Nonequilibrium_exchange::
Nonequilibrium_exchange(const string suffix, double A, double n, double E_a,
			int ivib, double alpha, double U, double D)
    : GeneralisedArrhenius(suffix, A, n, E_a),
      ivib_(ivib), alpha_(alpha), U_(U), Aj_(E_a), D_(D)
{
     model_ = "Nonequilibrium_exchange";
}

Nonequilibrium_exchange::
Nonequilibrium_exchange(ConfigParser &cfg,
			const string section, const string suffix)
    : GeneralisedArrhenius(cfg, section, suffix)
{
    Aj_ = E_a_;
    model_ = "Nonequilibrium_exchange";
    string alpha_s = "alpha" + suffix;
    if( ! cfg.parse_double(section, alpha_s, alpha_, 0.0) ) {
	cout << "Error reading in value " << alpha_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string ivib_s = "ivib" + suffix;
    

    if( ! cfg.parse_int(section, ivib_s, ivib_, -1) ) {
	cout << "Error reading in value for " << ivib_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( ivib_ < 0 ) {
	cout << "Error reading in value for ivib_f in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "The value should be listed as an integer, greater than or equal to 0.\n"
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string U_s = "U" + suffix;

    if( ! cfg.parse_double(section, U_s, U_, 0.0) ) {
	cout << "Error reading in value " << U_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string D_s = "D" + suffix;
    
    if( ! cfg.parse_double(section, D_s, D_, 0.0) ) {
	cout << "Error reading in value " << D_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
 
}

Nonequilibrium_exchange::
Nonequilibrium_exchange(const Nonequilibrium_exchange &n)
    : GeneralisedArrhenius(n.suffix_, n.A_, n.n_, n.E_a_),
      ivib_(n.ivib_), alpha_(n.alpha_), U_(n.U_), Aj_(n.Aj_),
      D_(n.D_)
{
     model_ = "Nonequilibrium_exchange";
}

Nonequilibrium_exchange::
~Nonequilibrium_exchange() {}


Nonequilibrium_exchange*
Nonequilibrium_exchange::
clone()
{
    return new Nonequilibrium_exchange(*this);
}


string
Nonequilibrium_exchange::
str() const
{
    ostringstream ost;
    ost << GeneralisedArrhenius::str()
	<< "alpha" << suffix_ << "                    = " << alpha_ << endl
	<< "ivib" << suffix_ << "                     = " << ivib_ << endl
	<< "U" << suffix_ << "                        = " << U_ << endl
	<< "D" << suffix_ << "                        = " << D_ << endl
	<< "A" << suffix_ << "                        = " << E_a_ << endl;
    return ost.str();
}

void
Nonequilibrium_exchange::
eval(const Gas_data &Q)
{
    // This call set k_
    GeneralisedArrhenius::eval(Q);
    double psi;
    if ( polyatomic_check(ivib_) ) {
        valarray<double> theta_vs;
        get_theta_vs(theta_vs,ivib_);
        psi = polyatomic_nonequilibrium_factor(Q, ivib_, theta_vs, alpha_, Aj_, U_, D_);
    } else
        psi = nonequilibrium_factor(Q, ivib_, get_theta_v(ivib_), alpha_, Aj_, U_, D_);
    // cout << "Nonequilibrium_exchange::eval() ivib= " << ivib_ << endl;
//     cout << "k_eq= " << k_ << endl;
    k_ = psi * k_;
    // cout << "psi= " << psi << " k= " << k_ << endl;
}

double nonequilibrium_factor(const Gas_data &Q, int ivib, double theta_v,
			     double alpha, double A, double U, double D)
{
    const double small_barrier = 1.0e-50;
    double Tvib = Q.T_vib[ivib];
    double T = Q.T;
    const double TOL = 1.0e-7;
    bool almost_eq = (fabs(Tvib - T) < TOL);
    if( almost_eq ) {
	return 1.0;
    }
    // 1. Calculate pseudo-temperatures
    double gamma_inv, gamma, T0, T_star;

   //  cout << "NONEQUILIBRIUM FACTOR: ivib= " << ivib << endl;
//     cout << "Tvib= " << Tvib << endl;
//     cout << "T= " << T << endl;

    // FIX WHEN U < 0.0
    if( U < 0.0 ) {
	gamma_inv = 1.0/Tvib - 1.0/T;
	gamma = 1.0/gamma_inv;
	T0 = Tvib;
	T_star = T;
    }
    else {
	gamma_inv = 1.0/Tvib - 1.0/T - 1.0/U;
	gamma = 1.0/gamma_inv;
	double T0_inv = 1.0/Tvib - 1.0/U;
	T0 = 1.0/T0_inv;
	double T_star_inv = 1.0/T - 1.0/U;
	T_star = 1.0/T_star_inv;
    }

   //  cout << "gamma= " << gamma << " T0= " << T0 << " T_star= " << T_star << endl;
//     cout << "alpha= " << alpha << " A= " << A << " U= " << U << " D= " << D << endl;
    
    double Q_D_T = vib_partition_function(theta_v, D, T);
    double Q_D_Tvib = vib_partition_function(theta_v, D, Tvib);
    double Q_alph_gam = vib_partition_function(theta_v, alpha*A, gamma);
    double Q_D_T0 = vib_partition_function(theta_v, D, T0);
    double Q_alph_T0 = vib_partition_function(theta_v, alpha*A, T0);
    double Q_alph_U;
    if( U < 0.0 ) {
	Q_alph_U = alpha*A/(theta_v*PC_k_SI);
    }
    else {
	Q_alph_U = vib_partition_function(theta_v, alpha*A, -U);
    }
    double Q_D_Ts = vib_partition_function(theta_v, D, T_star);
    if( A <= small_barrier ) {
	return ((Q_D_T*Q_D_T0)/(Q_D_Tvib*Q_D_Ts));
    }
    double Q_alph_Ts = vib_partition_function(theta_v, alpha*A, T_star);

    // cout << "Q_alph_gam= " << Q_alph_gam << " Q_D_T0= " << Q_D_T0 << " Q_alph_T0= " << Q_alph_T0 << endl;
//     cout << "Q_alph_U= " << Q_alph_U << " Q_D_Ts= " << Q_D_Ts << " Q_alph_Ts= " << Q_alph_Ts << endl;
//     cout << "exp= " << exp(-alpha*A/(PC_k_SI*T)) << endl;

    double numer = exp(-alpha*A/(PC_k_SI*T))*Q_alph_gam + Q_D_T0 - Q_alph_T0;
    double denom = exp(-alpha*A/(PC_k_SI*T))*Q_alph_U + Q_D_Ts - Q_alph_Ts;
    double psi = (Q_D_T/Q_D_Tvib)*(numer/denom);
    // Unfortunately at small values of temperature, the values for
    // Q_D_T0 and Q_alph_T0 are essentially the same (to machine precision).
    // Similarly, are Q_D_Ts and Q_alph_Ts essentially the same. Also, the
    // exponential term is very small (<1.0e-50).  All of these factors conspire
    // to give numer = 0.0 and denom = 0.0 which leads to psi = nan.
    // The fix is to set psi = 1.0 in these instances.  This will be reasonable
    // because this occure at low temperatures when T does not differ greatly from Tvib.
    // In these case the deviation form equiilibrium won't  be very strong and so
    // a value of psi = 1.0 will be ok.

    if( isnan(psi) || isinf(psi) )
	psi = 1.0;
    return psi;

}

double polyatomic_nonequilibrium_factor(const Gas_data &Q, int ivib,
                                        valarray<double> theta_vs,
                                        double alpha, double A, double U, double D)
{
    const double small_barrier = 1.0e-50;
    double Tvib = Q.T_vib[ivib];
    double T = Q.T;
    const double TOL = 1.0e-7;
    bool almost_eq = (fabs(Tvib - T) < TOL);
    if( almost_eq ) {
	return 1.0;
    }
    // 1. Calculate pseudo-temperatures
    double gamma_inv, gamma, T0, T_star;

   //  cout << "NONEQUILIBRIUM FACTOR: ivib= " << ivib << endl;
//     cout << "Tvib= " << Tvib << endl;
//     cout << "T= " << T << endl;

    // FIX WHEN U < 0.0
    if( U < 0.0 ) {
	gamma_inv = 1.0/Tvib - 1.0/T;
	gamma = 1.0/gamma_inv;
	T0 = Tvib;
	T_star = T;
    }
    else {
	gamma_inv = 1.0/Tvib - 1.0/T - 1.0/U;
	gamma = 1.0/gamma_inv;
	double T0_inv = 1.0/Tvib - 1.0/U;
	T0 = 1.0/T0_inv;
	double T_star_inv = 1.0/T - 1.0/U;
	T_star = 1.0/T_star_inv;
    }

   //  cout << "gamma= " << gamma << " T0= " << T0 << " T_star= " << T_star << endl;
//     cout << "alpha= " << alpha << " A= " << A << " U= " << U << " D= " << D << endl;
    
    double Q_D_T = polyatomic_vib_partition_function(theta_vs, D, T);
    double Q_D_Tvib = polyatomic_vib_partition_function(theta_vs, D, Tvib);
    double Q_alph_gam = polyatomic_vib_partition_function(theta_vs, alpha*A, gamma);
    double Q_D_T0 = polyatomic_vib_partition_function(theta_vs, D, T0);
    double Q_alph_T0 = polyatomic_vib_partition_function(theta_vs, alpha*A, T0);
    double Q_alph_U;
    if( U < 0.0 ) {
	Q_alph_U = alpha*A/(theta_vs.sum()*PC_k_SI);
    }
    else {
	Q_alph_U = polyatomic_vib_partition_function(theta_vs, alpha*A, -U);
    }
    double Q_D_Ts = polyatomic_vib_partition_function(theta_vs, D, T_star);
    if( A <= small_barrier ) {
	return ((Q_D_T*Q_D_T0)/(Q_D_Tvib*Q_D_Ts));
    }
    double Q_alph_Ts = polyatomic_vib_partition_function(theta_vs, alpha*A, T_star);

    // cout << "Q_alph_gam= " << Q_alph_gam << " Q_D_T0= " << Q_D_T0 << " Q_alph_T0= " << Q_alph_T0 << endl;
//     cout << "Q_alph_U= " << Q_alph_U << " Q_D_Ts= " << Q_D_Ts << " Q_alph_Ts= " << Q_alph_Ts << endl;
//     cout << "exp= " << exp(-alpha*A/(PC_k_SI*T)) << endl;

    double numer = exp(-alpha*A/(PC_k_SI*T))*Q_alph_gam + Q_D_T0 - Q_alph_T0;
    double denom = exp(-alpha*A/(PC_k_SI*T))*Q_alph_U + Q_D_Ts - Q_alph_Ts;
    double psi = (Q_D_T/Q_D_Tvib)*(numer/denom);
    // Unfortunately at small values of temperature, the values for
    // Q_D_T0 and Q_alph_T0 are essentially the same (to machine precision).
    // Similarly, are Q_D_Ts and Q_alph_Ts essentially the same. Also, the
    // exponential term is very small (<1.0e-50).  All of these factors conspire
    // to give numer = 0.0 and denom = 0.0 which leads to psi = nan.
    // The fix is to set psi = 1.0 in these instances.  This will be reasonable
    // because this occure at low temperatures when T does not differ greatly from Tvib.
    // In these case the deviation form equiilibrium won't  be very strong and so
    // a value of psi = 1.0 will be ok.

    if( isnan(psi) || isinf(psi) )
	psi = 1.0;
    return psi;

}

double vib_partition_function(double theta_v, double Y, double T)
{
    double numer = 1.0 - exp(-Y/(PC_k_SI*T));
    double denom = 1.0 - exp(-theta_v/T);
    return (numer/denom);
}

double polyatomic_vib_partition_function(valarray<double> theta_vs, double Y, double T)
{
    double numer, denom;
    double Q_vib = 1.0;

    /* The limiting temperature Y needs to be shared amongst all modes */
    Y /= double(theta_vs.size());
    
    for ( size_t i=0; i<theta_vs.size(); i++ ) {
        numer = 1.0 - exp(-Y/(PC_k_SI*T));
        denom = 1.0 - exp(-theta_vs[i]/T);
        Q_vib *= (numer/denom);
    }
    return Q_vib;
}

ParkNonequilibrium::
ParkNonequilibrium( const string suffix, double A, double n, double E_a, int ivib, double s )
    : RateCoeffModel( "ParkNonequilibrium", suffix),
      A_( A ), n_( n ), E_a_( E_a ), ivib_( ivib ), s_( s ) {}

ParkNonequilibrium::
ParkNonequilibrium( ConfigParser &cfg, const string section,
		      const string suffix )
    : RateCoeffModel( "ParkNonequilibrium", suffix)
{
    string k = "k" + suffix;
    vector<double> val;
    vector<double> notfound;
    notfound.push_back(-1.0);
    if( ! cfg.parse_vector_of_doubles( section, k, val, notfound ) ) {
	cerr << "ParkNonequilibrium::ParkNonequilibrium() --- "
	     <<  __FILE__ << (__LINE__ - 2) << endl
	     << "Error reading " << k << " in [" << section << "] of: "
	     << cfg.file_name
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    // Check the value makes sense.
    if( val.size() == 3 ) { // That's a good start
	A_ = val[0];
	n_ = val[1];
	E_a_ = val[2];
    }
    else {
	cerr << "ParkNonequilibrium::ParkNonequilibrium() --- "
	     <<  __FILE__ << (__LINE__ - 8) << endl
	     << "Error with number of values listed for ParkNonequilibrium\n"
	     << "coefficient parameters.\n"
	     << "Three values are exepected: A  n  E_a in expression A T^n exp(-E_a/(kT))\n"
	     << "A size of value vector " << val.size() << " has been found.\n"
	     << "Perhaps the key " << k << " has not been listed in [" << section << "] "
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string ivib_s = "ivib" + suffix;

    if( ! cfg.parse_int(section, ivib_s, ivib_, -1) ) {
	cout << "Error reading in value for " << ivib_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( ivib_ < 0 ) {
	cout << "Error reading in value for ivib_f in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "The value should be listed as an integer, greater than or equal to 0.\n"
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    string s_s = "s" + suffix;
    
    if( ! cfg.parse_double(section, s_s, s_, -1.0) ) {
	cout << "Error reading in value for " << s_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( s_ < 0.0 || s_ > 1.0 ) {
	cout << "Error reading in value for " << s_s << " in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "The value should be a double between 0.0 and 1.0.\n"
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    return;
}

ParkNonequilibrium::
ParkNonequilibrium( const ParkNonequilibrium &g )
    : RateCoeffModel( g.model_, g.suffix_ ),
      A_( g.A_ ), n_( g.n_ ), E_a_( g.E_a_ ), ivib_( g.ivib_ ), s_( g.s_ ) {}

ParkNonequilibrium::
~ParkNonequilibrium() {}

ParkNonequilibrium*
ParkNonequilibrium::
clone()
{
    return new ParkNonequilibrium(*this);
}


string
ParkNonequilibrium::
str() const
{
    ostringstream ost;
    ost << setprecision(6) << showpoint;
    ost << "k" << suffix_ << "                        = "
	<< A_ << "  "
	<< n_ << "  "
	<< E_a_ <<  endl
	<< "ivib" << suffix_ << "                     = " << ivib_ << endl
	<< "s" << suffix_ << "                        = " << s_ << endl;
    return ost.str();
}

void
ParkNonequilibrium::
eval( const Gas_data &Q )
{
    // Geometric average of T_tr and T_ve is used as rate controlling temperature
    double T_ga = pow(Q.T,s_) * pow(Q.T_vib[ivib_],(1-s_));
    k_ = A_ * pow(T_ga, n_) * exp(- E_a_ / (PC_k_SI * T_ga) );
}

RadiativeDecay::
RadiativeDecay( const string suffix, double lambda, double tau )
    : RateCoeffModel( "RadiativeDecay", suffix),
      lambda_( lambda ), tau_( tau ) {}

RadiativeDecay::
RadiativeDecay( ConfigParser &cfg, const string section,
		const string suffix )
    : RateCoeffModel( "RadiativeDecay", suffix)
{
    string k = "k" + suffix;
    vector<double> val;
    vector<double> notfound;
    notfound.push_back(-1.0);
    if( ! cfg.parse_vector_of_doubles( section, k, val, notfound ) ) {
	cerr << "RadiativeDecay::RadiativeDecay() --- "
	     <<  __FILE__ << (__LINE__ - 2) << endl
	     << "Error reading " << k << " in [" << section << "] of: "
	     << cfg.file_name
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    // Check the value makes sense.
    if( val.size() == 2 ) { // That's a good start
	lambda_ = val[0];
	tau_ = val[1];
    }
    else {
	cerr << "RadiativeDecay::RadiativeDecay() --- "
	     <<  __FILE__ << (__LINE__ - 8) << endl
	     << "Error with number of values listed for RadiativeDecay\n"
	     << "coefficient parameters.\n"
	     << "Three values are exepected: lambda tau in expression lambda/tau \n"
	     << "A size of value vector " << val.size() << " has been found.\n"
	     << "Perhaps the key " << k << " has not been listed in [" << section << "] "
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    return;
}

RadiativeDecay::
RadiativeDecay( const ParkNonequilibrium &g )
    : RateCoeffModel( g.model_, g.suffix_ ),
      lambda_( g.lambda_ ), tau_( g.tau_ ) {}

RadiativeDecay::
~RadiativeDecay() {}

RadiativeDecay*
RadiativeDecay::
clone()
{
    return new RadiativeDecay(*this);
}


string
RadiativeDecay::
str() const
{
    ostringstream ost;
    ost << setprecision(6) << showpoint;
    ost << "k" << suffix_ << "                        = "
	<< lambda_ << "  "
	<< tau_ << endl;
    return ost.str();
}

void
RadiativeDecay::
eval( const Gas_data &Q )
{
    // Simply transition probability (1/tau) times the escape factor (lambda)
    k_ = lambda_ / tau_;
}
