// Author: Brendan T. O'Flaherty
//  based on perfect-thermal-behaviour.cxx
// Date: 10-Aug-2009
// Place: Brisbane, Qld., A.

#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "../../nm/source/functor.hh"

#include "CEA-Cp-functor.hh"
#include "CEA-h-functor.hh"
#include "CEA-s-functor.hh"
#include "real-thermal-behaviour.hh"

using namespace std;

Real_thermal_behaviour::
Real_thermal_behaviour(lua_State *L)
    : Thermal_behaviour_model()
{
    T_COLD_ = get_positive_number(L, LUA_GLOBALSINDEX, "T_COLD");

    lua_getglobal(L, "species");

    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Real_thermal_behaviour::Real_thermal_behaviour():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    int nsp = lua_objlen(L, -1);

    for ( int isp = 0; isp < nsp; ++isp ) {
	vector<Univariate_functor*> Cp;
	vector<Univariate_functor*> h;
	vector<Univariate_functor*> s;
	vector<double> breaks;
	double T_low, T_high;

	lua_rawgeti(L, -1, isp+1);
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Real_thermal_behaviour::Real_thermal_behaviour():\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	double M = get_positive_value(L, -1, "M");
	M_.push_back(M);
	R_.push_back(PC_R_u/M);

	lua_getfield(L, -1, "CEA_coeffs");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Real_thermal_behaviour::Real_thermal_behaviour():\n";
	    ost << "Error locating 'CEA_coeffs' table for species: " << sp << endl;
	    input_error(ost);
	}
	for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i);
	    T_low = get_positive_number(L, -1, "T_low");
	    T_high = get_positive_number(L, -1, "T_high");

	    Cp.push_back(new CEA_Cp_functor(L, R_.back()));
	    h.push_back(new CEA_h_functor(L, (*Cp.back())(T_low), (*Cp.back())(T_high), R_.back()));
	    s.push_back(new CEA_s_functor(L, (*Cp.back())(T_low), (*Cp.back())(T_high), R_.back()));
	    breaks.push_back(T_low);
	    lua_pop(L, 1);
	}
	breaks.push_back(T_high);
	lua_pop(L, 2); // pop coeffs, species

	Cp_.push_back(new Segmented_functor(Cp, breaks));
	h_.push_back(new Segmented_functor(h, breaks));
	s_.push_back(new Segmented_functor(s, breaks));
	
	for ( size_t i = 0; i < Cp.size(); ++i ) {
	    delete Cp[i];
	    delete h[i];
	    delete s[i];
	} 
    }
}

Real_thermal_behaviour::
~Real_thermal_behaviour()
{
    for ( size_t isp = 0; isp < Cp_.size(); ++isp ) {
	delete Cp_[isp];
	delete h_[isp];
	delete s_[isp];
    }
}

int
Real_thermal_behaviour::
s_decode_conserved_energy(Gas_data &Q, const vector<double> &rhoe)
{
    return tbm_decode_conserved_energy(Q.e, rhoe, Q.rho);
}

int
Real_thermal_behaviour::
s_encode_conserved_energy(const Gas_data &Q, vector<double> &rhoe)
{
    return tbm_encode_conserved_energy(rhoe, Q.e, Q.rho);
}

double
Real_thermal_behaviour::
s_dhdT_const_p(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return tbm_dhdT_const_p(Cp_, Q.massf, Q.T);
}

int
Real_thermal_behaviour::
s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_)
{
    // given a correct value for e and rho
    // evaluate a new value for T.
    const int MAX_ATTEMPTS = 20;
    const int MAX_OSCIL = 3;
    const double TOL = 1.0e-6; // Experience shows this is a good value
    const bool use_T_guess = true;
    
    double e_given = Q.e[0];
    double T_0 = 0.0;
    double T_1 = 0.0;
    double T_prev = 0.0;
    double alpha = 0.0;

    // Because the CEA curve fits cut off at 200.0, the iteration
    // doesn't always work well when starting below this value.
    // So when that's the case, we use our crude starting guess.
    if ( use_T_guess && Q.T[0] > 200.0 )
	T_0 = Q.T[0];
    else
	T_0 = Q.e[0] / 717.0; // C_v for "normal" air as a crude guess
    
    // We use a Newton solver to iterate for temperature
    // 1. Evalute f_0 based on guess
    EOS_->eval_pressure(Q);
    double f_0 = zero_function(Q, EOS_, e_given, T_0);
    double dfdT = deriv_function(Q, EOS_, T_0);

    if( fabs( f_0 ) < TOL )
	return SUCCESS;

    T_prev = T_0;
    
    // DFP: Seem to encounter successive oscillations around T = 6000 K
    //      Introducing a flag to count them and break the iterative
    //      loop after a set number oscillations (currently three).
    
    int oscil_count = 0;
    int attempt;
    for ( attempt = 0; attempt < MAX_ATTEMPTS; ++attempt ) {
	if ( fabs(dfdT) < 1.0e-12 ) {
	    cout << "WARNING - Real_thermal_behaviour::solve_for_T()\n";
	    cout << "Nearly zero derivative, dfdT= " << dfdT << endl;
	    return FAILURE;
	}
	T_1 = T_0 - f_0/dfdT;
	Q.T[0] = T_1;
	EOS_->eval_pressure(Q);
	if ( T_1 < T_COLD_ ) {
	    Q.T[0] = T_COLD_;
	    s_eval_energy(Q, EOS_);
	    return SUCCESS;
	}
	if ( fabs(T_0 - T_1) < TOL ) {
	    Q.T[0] = T_1;
	    Q.e[0] = e_given;
	    return SUCCESS;
	}
	if ( fabs(T_1 - T_prev) < 0.01*TOL ) {
	    // This is when we hit a symmetric oscillation
            if ( ++oscil_count == MAX_OSCIL ) break;
	    alpha = fabs((double) rand() / (double) RAND_MAX);
	    T_prev = T_0;
	    T_1 = T_0 - alpha*f_0/dfdT;
	    T_0 = T_1;
	    f_0 = zero_function(Q, EOS_, e_given, T_0);
	    dfdT = deriv_function(Q, EOS_, T_0);
	}
	else {
	    T_prev = T_0;
	    T_0 = T_1;
	    f_0 = zero_function(Q, EOS_, e_given, T_0);
	    dfdT = deriv_function(Q, EOS_, T_0);
	}
    }

    if( test_T_for_polynomial_breaks(T_1) ) {
	Q.T[0] = T_1;
	Q.e[0] = e_given;
    }
    
    // If we get this far, then the iterations did not converge.
    // We should print a warning and return -1.0 as temperature.
    // DFP: Customised the error message for the case where 
    //      successive symmetric oscillations has been determined.
    cout << "WARNING - ThermallyRealGasMix::solve_for_T() \n";
    if ( oscil_count == MAX_OSCIL ) {
        cout << oscil_count << " symmetric oscillations were encountered by the Newton solver.\n";
	cout << "In addition, the temperature does not correspond to a CEA polynomial break.\n";
    }
    else {
	cout << "The Newton solver did not converge after " << attempt << " iterations.\n";
    }
    cout << "T_1= " << T_1 << " T_0= " << T_0 << " f_0= " << f_0 << " dfdT= " << dfdT << endl;
    if ( use_T_guess ) {
	cout << "T_guess - from gas_data = " << Q.T[0] << endl;
    }
    else {
        cout << "T_guess - from e[0], Cv = " <<  Q.e[0] / 717.0 << endl;;
    }
    cout << "Gas state...\n";
    Q.print_values();
    cout << "Bailing out!\n";
    exit(ITERATION_ERROR);
}

double
Real_thermal_behaviour::
s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    // Reference:
    // Cengel and Boles (2002)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // 4th Ed.
    // McGraw Hill
    // Equation 11-46 on p. 617
    
    double nu = 1.0/Q.rho;
    double dnudT = -(nu*nu)/EOS_->dTdrho_const_p(Q, status);
    double dpdnu = -1.0/(nu*nu)*EOS_->dpdrho_const_T(Q, status);
    double Cv = 0.0;
    for ( size_t isp = 0; isp < Cp_.size(); ++isp ) {
	Cv += Q.massf[isp]*((*Cp_[isp])(Q.T[0]));
    }
    Cv += Q.T[0]*dnudT*dnudT*dpdnu;
    return Cv;
}

int
Real_thermal_behaviour::
s_eval_energy(Gas_data &Q, Equation_of_state *EOS_)
{
    double e = 0.0;
    for ( size_t isp = 0; isp < Cp_.size(); ++isp ) {
	e += Q.massf[isp]*s_eval_energy_isp(Q, EOS_, isp);
    }
    Q.e[0] = e;
    return SUCCESS;
}

double
Real_thermal_behaviour::
s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    // Reference:
    // Cengel and Boles (2002)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // 4th Ed.
    // McGraw Hill
    // Equation 11-29 on p. 614

    double h = s_eval_enthalpy_isp(Q, EOS_, isp);
    double pv = EOS_->prho_ratio(Q, isp);
    return h - pv;
}

double
Real_thermal_behaviour::
s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    // Reference:
    // Cengel and Boles (2002)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // 4th Ed.
    // McGraw Hill
    // Equation 11-35 on p. 615

    int status;
    double nu = 1.0/Q.rho;
    double dnudT = -(nu*nu)/EOS_->dTdrho_const_p(Q, status);
    double M = calculate_molecular_weight(Q.massf, M_);
    return (*h_[isp])(Q.T[0]) + (nu - Q.T[0]*dnudT)*Q.p*(M/M_[isp]);
}

double
Real_thermal_behaviour::
s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    // Reference:
    // Cengel and Boles (2002)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // 4th Ed.
    // McGraw Hill
    // Equation 11-40 on p. 615

    int status;
    double nu = 1.0/Q.rho;
    double dnudT = -(nu*nu)/EOS_->dTdrho_const_p(Q, status);
    double M = calculate_molecular_weight(Q.massf, M_);

    return (*s_[isp])(Q.T[0]) - dnudT*Q.p*(M/M_[isp]);
}

double
Real_thermal_behaviour::
zero_function(Gas_data &Q, Equation_of_state *EOS_, double e_given, double T)
{
    // given correct energy
    // evaluate the zero function using temperature
    // f(T) = e_given - e(T) = 0
    Q.T[0] = T;
    s_eval_energy(Q, EOS_);
    return Q.e[0] - e_given;
}

double
Real_thermal_behaviour::
deriv_function(Gas_data &Q, Equation_of_state *EOS_, double T)
{
    // return dedT_const_v
    Q.T[0] = T;
    int status;
    return dedT_const_v(Q, EOS_, status);
}

bool
Real_thermal_behaviour::
test_T_for_polynomial_breaks(double T)
{
    // Sometimes the iterations run into trouble because of badly formed
    // polynomials near the breaks.
    // It turns out the the CEA polynomials are not that continuous on the
    // the finer detail.

    if( T >= 999.99 && T <= 1000.01 )
	return true;
    if( T >= 5999.99 && T <= 6000.01 )
	return true;

    return false;
}
