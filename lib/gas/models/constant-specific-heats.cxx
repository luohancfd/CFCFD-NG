// Author: Rowan J. Gollan
// Version: 24-May-2008
//            Initial coding
//

#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "constant-specific-heats.hh"
#include "physical_constants.hh"

using namespace std;

Constant_specific_heats::
Constant_specific_heats(lua_State *L)
    : Thermal_behaviour_model()
{
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Constant_specific_heats::Constant_specific_heats():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    int nsp = lua_objlen(L, -1);
    
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1); // A Lua list is offset one from the C++ vector index
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring the specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Constant_specific_heats::Constant_specific_heats():\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	double M = get_positive_value(L, -1, "M");
	double R = PC_R_u/M;
	double gamma = get_positive_value(L, -1, "gamma");
	double Cv = R / (gamma - 1.0);
	Cv_.push_back(Cv);
	double Cp = R + Cv;
	Cp_.push_back(Cp);
	double e_zero = get_value(L, -1, "e_zero");
	e_zero_.push_back(e_zero);
	double q = get_value(L, -1, "q");
	q_.push_back(q);
	s_ref_.push_back(Cp*log(PC_T_ref));
	
	lua_pop(L, 1); // pop "sp" off stack
    }
    lua_pop(L, 1); // pop "species" off stack
}

Constant_specific_heats::
~Constant_specific_heats() {}

double
Constant_specific_heats::
s_dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    status = SUCCESS;
    return mass_average(Q.massf, Cp_);
}

double
Constant_specific_heats::
s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    status = SUCCESS;
    return mass_average(Q.massf, Cv_);
}

int
Constant_specific_heats::
s_eval_energy(Gas_data &Q, Equation_of_state *EOS_)
{
    double Cv = mass_average(Q.massf, Cv_);
    double e_zero = mass_average(Q.massf, e_zero_);
    double q = mass_average(Q.massf, q_);
    Q.e[0] = Cv*Q.T[0] - e_zero - q;
    return SUCCESS;
}

int
Constant_specific_heats::
s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_)
{
    double Cv = mass_average(Q.massf, Cv_);
    double e_zero = mass_average(Q.massf, e_zero_);
    double q = mass_average(Q.massf, q_);
    Q.T[0] = (Q.e[0] + e_zero + q)/Cv;
    return SUCCESS;
}

double
Constant_specific_heats::
s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    return Cv_[isp]*Q.T[0];
}

double
Constant_specific_heats::
s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    double R = Cp_[isp] - Cv_[isp];
    return Cp_[isp]*log(Q.T[0]/PC_T_ref) - R*log(Q.p/PC_P_atm) + s_ref_[isp];
}

double
Constant_specific_heats::
s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    return Cp_[isp]*Q.T[0] - e_zero_[isp] - q_[isp];
}
