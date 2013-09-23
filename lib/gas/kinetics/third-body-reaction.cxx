// Author: Rowan J. Gollan
// Date: 29-Mar-2009
// Place: Poquoson, Virginia, USA

#include <iostream>
#include "../../util/source/lua_service.hh"
#include "third-body-reaction.hh"

using namespace std;

Third_body_reaction::
Third_body_reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower)
    : Normal_reaction(L, g, T_upper, T_lower), conc_(0.0),
      conc_just_computed_(false)
{
    read_table_as_map(L, -1, "efficiencies", efficiencies_);
}

Third_body_reaction::
~Third_body_reaction() {}

double
Third_body_reaction::
s_compute_forward_rate(const valarray<double> &y)
{
    double wfN = Normal_reaction::s_compute_forward_rate(y);
    compute_third_body_concentration(y);
    return conc_*wfN;
}

double
Third_body_reaction::
s_compute_backward_rate(const valarray<double> &y)
{
    double wbN = Normal_reaction::s_compute_backward_rate(y);
    if ( !conc_just_computed_ )
	compute_third_body_concentration(y);
    conc_just_computed_ = false;
    return conc_*wbN;
}

void
Third_body_reaction::
compute_third_body_concentration(const valarray<double> &y)
{
    conc_ = 0.0;
    map<int, double>::const_iterator it;
    for ( it = efficiencies_.begin(); it != efficiencies_.end(); ++it ) {
	conc_ += it->second * y[it->first];
    }
    conc_just_computed_ = true;
}

Reaction* create_Third_body_reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower)
{
    return new Third_body_reaction(L, g, T_upper, T_lower);
}
