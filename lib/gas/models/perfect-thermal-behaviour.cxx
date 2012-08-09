// Author: Rowan J. Gollan
// Date: 03-Nov-2008
// Place: Hampton, Virginia, USA

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
#include "perfect-thermal-behaviour.hh"

using namespace std;

Perfect_thermal_behaviour::
Perfect_thermal_behaviour(lua_State *L)
    : Thermal_behaviour_model()
{
    T_COLD_ = get_positive_number(L, LUA_GLOBALSINDEX, "T_COLD");

    lua_getglobal(L, "species");

    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Perfect_thermal_behaviour::Perfect_thermal_behaviour():\n";
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
	    ost << "Perfect_thermal_behaviour::Perfect_thermal_behaviour():\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	double M = get_positive_value(L, -1, "M");
	R_.push_back(PC_R_u/M);

	lua_getfield(L, -1, "CEA_coeffs");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Perfect_thermal_behaviour::Perfect_thermal_behaviour():\n";
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

Perfect_thermal_behaviour::
~Perfect_thermal_behaviour()
{
    for ( size_t isp = 0; isp < Cp_.size(); ++isp ) {
	delete Cp_[isp];
	delete h_[isp];
	delete s_[isp];
    }
}

int
Perfect_thermal_behaviour::
s_decode_conserved_energy(Gas_data &Q, const vector<double> &rhoe)
{
    return tbm_decode_conserved_energy(Q.e, rhoe, Q.rho);
}

int
Perfect_thermal_behaviour::
s_encode_conserved_energy(const Gas_data &Q, vector<double> &rhoe)
{
    return tbm_encode_conserved_energy(rhoe, Q.e, Q.rho);
}

double
Perfect_thermal_behaviour::
s_dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    status = SUCCESS;
    return tbm_dhdT_const_p(Cp_, Q.massf, Q.T);
}

double
Perfect_thermal_behaviour::
s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    double Cv = 0.0;
    for ( size_t isp = 0; isp < Cp_.size(); ++isp ) {
	Cv += Q.massf[isp]*((*Cp_[isp])(Q.T[0]) - R_[isp]);
    }
    status = SUCCESS;
    return Cv;
}

int
Perfect_thermal_behaviour::
s_eval_energy(Gas_data &Q, Equation_of_state *EOS_)
{
    int status = SUCCESS;
    double e = 0.0;
    for ( size_t isp = 0; isp < Cp_.size(); ++isp ) {
	e += Q.massf[isp]*s_eval_energy_isp(Q, EOS_, isp);
    }
    Q.e[0] = e;
    return status;
}

int
Perfect_thermal_behaviour::
s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_)
{
    const int MAX_ATTEMPTS = 30;
    const int MAX_OSCIL = 3;
    const double TOL = 1.0e-6; // Experience shows this is a good value
    const bool use_T_guess = true;
    const double rel_factor = 0.1; // For use when iterations run into trouble

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

    double T_guess = T_0;
    
    // We use a Newton solver to iterate for temperature
    // 1. Evalute f_0 based on guess
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
	    cout << "WARNING - Perfect_thermal_behaviour::solve_for_T()\n";
	    cout << "Nearly zero derivative, dfdT= " << dfdT << endl;
	    return FAILURE;
	}
	T_1 = T_0 - f_0/dfdT;
	
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
	return SUCCESS;
    }

    // If we get this far, then the iterations did not converge.
    // One more thing to try.
    // Continue iterating but only use a fraction of the suggested
    // update.
    
    for ( attempt = 0; attempt < MAX_ATTEMPTS; ++attempt ) {
	if ( fabs(dfdT) < 1.0e-12 ) {
	    cout << "WARNING - Perfect_thermal_behaviour::solve_for_T()\n";
	    cout << "Nearly zero derivative, dfdT= " << dfdT << endl;
	    return FAILURE;
	}
	T_1 = T_0 - rel_factor*f_0/dfdT;
	
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

    // We should print a warning and return -1.0 as temperature.
    // DFP: Customised the error message for the case where 
    //      successive symmetric oscillations has been determined.
    cout << "WARNING - Perfect_thermal_behaviour::solve_for_T() \n";
    if ( oscil_count == MAX_OSCIL ) {
        cout << oscil_count << " symmetric oscillations were encountered by the Newton solver.\n";
	cout << "In addition, the temperature does not correspond to a CEA polynomial break.\n";
    }
    else {
	cout << "The Newton solver did not converge after " << attempt << " iterations.\n";
    }
    cout << "T_1= " << T_1 << " T_0= " << T_0 << " f_0= " << f_0 << " dfdT= " << dfdT << endl;
    if ( use_T_guess ) {
	cout << "T_guess - from gas_data = " << T_guess << endl;
    }
    else {
        cout << "T_guess - from e[0], Cv = " <<  T_guess << endl;;
    }
    cout << "Gas state...\n";
    Q.e[0] = e_given;
    Q.print_values();
    cout << "Bailing out!\n";
    exit(ITERATION_ERROR);
}

double
Perfect_thermal_behaviour::
zero_function(Gas_data &Q, Equation_of_state *EOS_, double e_given, double T)
{
    Q.T[0] = T;
    s_eval_energy(Q, EOS_);
    return Q.e[0] - e_given;
}

double
Perfect_thermal_behaviour::
deriv_function(Gas_data &Q, Equation_of_state *EOS_, double T)
{
    Q.T[0] = T;
    int status;
    return dedT_const_v(Q, EOS_, status);
}


bool
Perfect_thermal_behaviour::
test_T_for_polynomial_breaks(double T)
{
    // Sometimes the iterations run into trouble because of badly formed
    // polynomials near the breaks.
    // It turns out the the CEA polynomials are not that continuous on the
    // the finer detail.
    // This error should only occur if all break points are the same, and 
    // we have landed on top of one.

    double dT=0.01;
    for (size_t i=0; i<h_[0]->size(); ++i) {
	if (T >= h_[0]->get_break(i)-dT && T <= h_[0]->get_break(i)+dT) {
	    return true;
	}
    }

    return false;
}

double
Perfect_thermal_behaviour::
s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    double h = (*h_[isp])(Q.T[0]);
    double e = h - EOS_->prho_ratio(Q, isp);

    return e;
}

double
Perfect_thermal_behaviour::
s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    return (*h_[isp])(Q.T[0]);
}

double
Perfect_thermal_behaviour::
s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    return (*s_[isp])(Q.T[0]) - R_[isp]*log(Q.p/PC_P_atm);
}
