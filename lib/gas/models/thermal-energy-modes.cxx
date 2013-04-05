// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#include <set>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"
#include "thermal-energy-modes.hh"
#include "gas-model.hh"

using namespace std;

Thermal_energy_mode::
Thermal_energy_mode( string name, vector<Chemical_species*> &species, lua_State *L )
: name_( name)
{
    // Governing temperature
    iT_ = get_int( L, -1, "iT" );
    
    // Thermal components list - allows gathing of components from species
    lua_getfield( L, -1, "components" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Thermal_energy_mode::Thermal_energy_mode()\n";
	ost << "Error locating 'components' table" << endl;
	input_error(ost);
    }
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	component_names_.push_back( luaL_checkstring(L, -1) );
	lua_pop(L, 1 );
    }
    
    lua_pop(L, 1);	// pop 'components'

    // Check validity of component names then create components
    cout << "- Creating a Thermal_energy_mode '" << name_ << "' with components:" << endl;
    string valid_mode_types[] = { "translation", "rotation", "vibration", "electronic", "internal" };
    for ( size_t ic=0; ic<component_names_.size(); ++ic ) {
    	bool valid_type = false;
    	string name = component_names_[ic];
    	cout << "  * " << name << " - { ";
	for ( int itype=0; itype<5; ++itype ) {
	    if ( name.find( valid_mode_types[itype] )!=string::npos ) {
		valid_type = true;
		break;
	    }
	    // else -> continue
	}
	if ( !valid_type ) {
	    ostringstream ost;
	    ost << "Thermal_energy_mode::Thermal_energy_mode():\n";
	    ost << "component name: " << name <<" is not valid.\n";
	    input_error(ost);
	}
	// else -> use modes types to test for inclusion as components
	for ( size_t isp=0; isp<species.size(); ++isp ) {
	    Chemical_species * X = species[isp];
	    // determine if this species is a heavy particle or not
	    string hp = "XX";
	    if ( X->get_name()!="e_minus" ) hp = "hp";
	    for ( int iem=0; iem<X->get_n_modes(); ++iem ) {
		Species_energy_mode * M = X->get_mode_pointer(iem);
		// Test for: (1) generic types eg "all-vibration" or "hp-translation"
		//           (2) specifically this species plus type eg "N2-vibration"
		if ( name.find( "all-"+M->get_type() )!=string::npos || 
		     name.find( hp + "-" + M->get_type() )!=string::npos ||
		     name.find( X->get_name() + "-" + M->get_type() )!=string::npos ) {
		    M->set_iT( iT_ );
		    components_.push_back( M );
		    // FIXME: Remove after testing
		    cout << X->get_name() << "-" << M->get_type() << ", ";
		    sp_idx_.push_back(M->get_isp());
		}
	    }
	}
	cout << " } " << endl;
    }
    // Convert sp_idx vector to a set to remove duplicate entries
    set<int> sp_idx2(sp_idx_.begin(), sp_idx_.end());
    // Assign set back to vector
    sp_idx_.assign(sp_idx2.begin(), sp_idx2.end());
}

Thermal_energy_mode::~Thermal_energy_mode()
{
    // components are just pointers - no deletes necessary for now
}

double
Thermal_energy_mode::
mode_massf(const Gas_data &Q)
{
    double sum = 0.0;
    for ( size_t i = 0; i < sp_idx_.size(); ++i ) {
	int isp = sp_idx_[i];
	sum += Q.massf[isp];
    }
    return sum;
}

double
Thermal_energy_mode::
s_decode_conserved_energy(Gas_data &Q, double rhoe)
{
    double modef = mode_massf(Q);
    if ( modef > DEFAULT_MIN_MASS_FRACTION ) {
	return rhoe/(Q.rho*modef);
    }
    else {
	return rhoe/Q.rho;
    }
}

double
Thermal_energy_mode::
s_encode_conserved_energy(const Gas_data &Q)
{
    return Q.rho*mode_massf(Q)*Q.e[iT_];
}

double
Thermal_energy_mode::
s_eval_dedT(Gas_data &Q)
{
    double dedT = 0.0;
    double modef = mode_massf(Q);
    for ( size_t ic=0; ic<components_.size(); ++ic ) {
	double Cv_i = components_[ic]->eval_Cv(Q);
	int isp = components_[ic]->get_isp();
	dedT += ( modef > DEFAULT_MIN_MASS_FRACTION ) ? (Q.massf[isp]/modef)*Cv_i : Cv_i;
    }
    
    return dedT;
}

double
Thermal_energy_mode::
s_eval_energy(const Gas_data &Q)
{
    double e = 0.0;
    double modef = mode_massf(Q);
    for ( size_t ic=0; ic<components_.size(); ++ic ) {
	double e_c = components_[ic]->eval_energy(Q);
	int isp = components_[ic]->get_isp();
	e += ( modef > DEFAULT_MIN_MASS_FRACTION ) ? (Q.massf[isp]/modef)*e_c : e_c; 
    }
    return e;
}

Constant_Cv_energy_mode::
Constant_Cv_energy_mode( string name, vector<Chemical_species*> &species, lua_State *L  )
 : Thermal_energy_mode( name, species, L ) {}
 
double
Constant_Cv_energy_mode::
s_eval_temperature(Gas_data &Q)
{
    // Direct solve for T as e = Cv.T
    
    return Q.e[iT_]/s_eval_dedT(Q);
}

void
Constant_Cv_energy_mode::
s_test_derivatives( Gas_data &Q )
{
    cout << "Constant_Cv_energy_mode::s_test_derivatives()" << endl
         << "The derivative is constant..." << endl
         << "Cv = " << s_eval_dedT(Q) << endl;
         
    return;
}

Variable_Cv_energy_mode::
Variable_Cv_energy_mode( string name, vector<Chemical_species*> &species, lua_State *L  )
 : Thermal_energy_mode( name, species, L ) {
    // Temperature range
    T_min_ = get_positive_number( L, -1, "T_min" );
    T_max_ = get_positive_number( L, -1, "T_max" );
    
    // Iterative method
    string iterative_method = get_string( L, -1, "iterative_method" );
    if ( iterative_method!="NewtonRaphson") {
    	ostringstream oss;
    	oss << "Thermal_energy_mode::Thermal_energy_mode()" << endl
    	    << "The requested iterative method: " << iterative_method 
    	    << " is not available." << endl;
    	input_error( oss );
    }
    
    // Maximum iterations
    max_iterations_ = get_int( L, -1, "max_iterations" );
    
    // Convergence tolerance 
    convergence_tolerance_ = get_positive_number( L, -1, "convergence_tolerance" );
}

bool
Variable_Cv_energy_mode::
check_T_range( double T )
{
    return ( T > T_min_ && T < T_max_ ) ? true : false;
}

void
Variable_Cv_energy_mode::
impose_T_limits( double &T )
{
    if ( T < T_min_ ) {
    	T = T_min_;
    }
    else if ( T > T_max_ ) {
    	T = T_max_;
    }
}

bool
Variable_Cv_energy_mode::
check_convergence( double e, double e_given )
{
    double error = fabs( e - e_given ) / e_given;
    
    return ( error < convergence_tolerance_ ) ? true : false;
}

double
Variable_Cv_energy_mode::
f_zero( double e, double e_given )
{
    return e - e_given;
}

double
Variable_Cv_energy_mode::
s_eval_temperature(Gas_data &Q)
{
    // Newton-Raphson iterations to solve for Q.T[iT_]
    
    // 0. Save initial params
    double e_given = Q.e[iT_];
    double T_given = Q.T[iT_];
    
    // 1.  Obtain an intial guess for T
    double T_im1;
    if ( check_T_range(T_given) ) {
    	// old T is valid
    	T_im1 = T_given;
    }
    else {
    	// make a crude guess using the bisection method
        T_im1 = s_eval_temperature_bisection(Q,1.0);
    }

    // 2. Set e_im1_
    Q.T[iT_] = T_im1;
    double e_im1 = s_eval_energy( Q );

    // 3. Perform iterations
    double T_i = T_given, e_i;
    for ( int i=0; i<max_iterations_; ++i ) {
    	// 3.a New T guess
    	T_i = T_im1 - f_zero( e_im1, e_given ) / dfdT( Q );
    	if ( T_i < T_min_ ) {
    	    T_i = T_min_;
    	    if ( T_im1 == T_max_ ) {
//    	    	cout << "Variable_Cv_energy_mode::s_eval_temperature()" << endl
//    	    	     << "Temperature is oscillating between T_min: " << T_min_
//    	    	     << " and T_max: " << T_max_ << endl;
//                Q.print_values(false);
//                Q.T[iT_] = T_min_;
//                double e_min = s_eval_energy( Q );
//                Q.T[iT_] = T_max_;
//                double e_max = s_eval_energy( Q );
//                cout << "e_given = " << e_given << ", e(T_min) = " << e_min << ", e(T_max) = " << e_max << endl;
                break;
    	    }
    	}
    	else if ( T_i > T_max_ ) {
    	    T_i = T_max_;
    	    if ( T_im1 == T_min_ ) {
//    	    	cout << "Variable_Cv_energy_mode::s_eval_temperature()" << endl
//    	    	     << "Temperature is oscillating between T_min: " << T_min_
//    	    	     << " and T_max: " << T_max_ << endl;
//    	        Q.print_values(false);
//                Q.T[iT_] = T_min_;
//                double e_min = s_eval_energy( Q );
//                Q.T[iT_] = T_max_;
//                double e_max = s_eval_energy( Q );
//                cout << "e_given = " << e_given << ", e(T_min) = " << e_min << ", e(T_max) = " << e_max << endl;
                break;
    	    }
    	}
    	Q.T[iT_] = T_i;
    	// 3.b New e level
    	e_i = s_eval_energy( Q );
        //cout << "i = " << i << ", T_im1 = " << T_im1 << ", T_i = " << T_i << ", e_i = " << e_i << ", e_given = " << e_given << endl;
    	// 3.c Test for convergence
    	if ( check_convergence( e_i, e_given ) ) break;
    	// else -> prepare for next iteration
    	T_im1 = T_i; e_im1 = e_i;
    }
    
    if ( !check_T_range( T_i ) ) {
        // A symmetric oscillation was probably encountered
        T_i = s_eval_temperature_bisection(Q,1.0e-6);
    }

    return T_i;
}

double
Variable_Cv_energy_mode::
s_eval_temperature_bisection(Gas_data &Q, double tol)
{
    // Bisection method to solve for Q.T[iT_]
    double e_given = Q.e[iT_];
    Q.T[iT_] = T_min_;
    double e_min = s_eval_energy( Q );
    Q.T[iT_] = T_max_;
    double e_max = s_eval_energy( Q );

    if ( e_given >= e_max ) {
        cout << "Variable_Cv_energy_mode::s_eval_temperature_bisection()" << endl
             << "Maximum temperature limit exceeded!" << endl;
        return T_max_;
    }
    else if ( e_given <= e_min ) {
        cout << "Variable_Cv_energy_mode::s_eval_temperature_bisection()" << endl
             << "Minimum temperature limit exceeded!" << endl;
        return T_min_;
    }

    double T_left = T_min_;
    double T_right = T_max_;

    double T_mid=(T_left+T_right)/2.0;

    for(T_mid=(T_left+T_right)/2.0; fabs(T_left-T_mid) > tol; T_mid=(T_left+T_right)/2.0) {
        // cout << "T_left = " << T_left << ", T_right = " << T_right << endl;
        Q.T[iT_] = T_left;
        double f_left = s_eval_energy( Q ) - e_given;
        Q.T[iT_] = T_right;
        double f_right = s_eval_energy( Q ) - e_given;
        if (f_left*f_right <= 0.0) {
            T_right = T_mid; // use left interval
        } else {
            T_left = T_mid; // use right interval
        }
    }

    return T_mid;
}


void
Variable_Cv_energy_mode::
s_test_derivatives( Gas_data &Q )
{
    // Write f_zero and analytical and numerical evaluations of dfdT over 
    // the given temperature range to file
    
    // Save initial params
    double e_given = Q.e[iT_];
    double T_given = Q.T[iT_];
    
    ofstream fz_file;
    fz_file.open("f_zero.txt",ios::out);
    fz_file << setprecision(12) << scientific << showpoint;
    fz_file << "# Column 1: T (K)" << endl
            << "# Column 2: f_zero: (J/kg)" << endl
            << "# Column 3: dfdT-analytic (J/kg-K)" << endl
            << "# Column 4: dfdT-numeric (J/kg-K)" << endl;
    double T = T_min_;
    double dT = ( T_max_ - T_min_ ) / 1000.0;
    while ( T <= T_max_ ) {
    	Q.T[iT_] = T;
    	double f = f_zero( s_eval_energy( Q ), e_given );
    	double dfdT_a = dfdT(Q);
    	Q.T[iT_] = T*0.99999;
    	double f_im1 = f_zero( s_eval_energy( Q ), e_given );
    	Q.T[iT_] = T*1.00001;
    	double f_ip1 = f_zero( s_eval_energy( Q ), e_given );
    	double dfdT_n = ( f_ip1 - f_im1 ) / (0.00002*T);
    	fz_file << setw(20) << T << setw(20) << f << setw(20) << dfdT_a << setw(20) << dfdT_n << endl;
    	T += dT;
    }
    fz_file.close();
    
    // Restore initial params
    Q.T[iT_] = T_given;
    Q.e[iT_] = e_given;
    
    return;
}

