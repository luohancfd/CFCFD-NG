// Author: Rowan J. Gollan
// Version: 17-Oct-2008
// Place: Hampton, Virginia, USA

#include <cmath>

#include "../../util/source/lua_service.hh"
#include "normal-reaction.hh"

using namespace std;

Normal_reaction::
Normal_reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower)
    : Reaction(L, g, T_upper, T_lower)
{
    read_table_as_map(L, -1, "f_coeffs", f_coeffs_);
    read_table_as_map(L, -1, "b_coeffs", b_coeffs_);
}

Normal_reaction::
~Normal_reaction() {}

double
Normal_reaction::
s_compute_forward_rate(const vector<double> &y)
{
    double val = 1.0;
    map<int, int>::const_iterator it;
    for ( it = f_coeffs_.begin(); it != f_coeffs_.end(); ++it ) {
	val *= pow(y[it->first], it->second);
    }
    return k_f()*val;
}

double
Normal_reaction::
s_compute_backward_rate(const vector<double> &y)
{
    double val = 1.0;
    map<int, int>::const_iterator it;
    for ( it = b_coeffs_.begin(); it != b_coeffs_.end(); ++it ) {
	val *= pow(y[it->first], it->second);
    }
    return k_b()*val;
}

Reaction* create_Normal_reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower)
{
    return new Normal_reaction(L, g, T_upper, T_lower);
}
