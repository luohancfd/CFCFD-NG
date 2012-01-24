// Author: Daniel F Potter
// Date: 07-Dec-2009

#include <cmath>
#include <sstream>
#include <iostream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "Park-nonequilibrium.hh"
#include "../models/physical_constants.hh"
#include "../models/chemical-species-library.hh"

using namespace std;

Park_nonequilibrium::
Park_nonequilibrium(lua_State *L, Gas_model &g)
    : Generalised_Arrhenius(L,g)
{
    // 1. Primary species data
    string p_name = get_string(L,-1,"p_name");
    string p_mode = get_string(L,-1,"p_mode");
    s_p_ = get_number(L, -1, "s_p");
    if ( s_p_ < 0.0 || s_p_ > 1.0 ) {
	ostringstream ost;
	ost << "Park_nonequilibrium::Park_nonequilibrium()" << endl
	    << "s_p was given as: " << s_p_ << " but is required to be between 0 and 1" << endl;
	input_error( ost );
    }
    iTp_ = get_library_species_pointer_from_name( p_name )->get_mode_pointer_from_type( p_mode )->get_iT();
    
    // 2. Secondary species data
    string q_name = get_string(L,-1,"q_name");
    string q_mode = get_string(L,-1,"q_mode");
    if ( q_name=="NA" ) {
    	if ( s_p_ != 1.0 ) {
    	    ostringstream ost;
    	    ost << "Park_nonequilibrium::Park_nonequilibrium()" << endl
    	        << "'q_name' and 'q_mode' are required when s_p is not 1" << endl;
    	    input_error( ost );
    	}
    	// just set iTq_ to iTp_ as a dummy value (will not be used)
    	iTq_ = iTp_;
    }
    else {
    	iTq_ = get_library_species_pointer_from_name( q_name )->get_mode_pointer_from_type( q_mode )->get_iT();
    }
    
}

Park_nonequilibrium::
Park_nonequilibrium(double A, double n, double E_a, double s_p, int iTp, int iTq)
    : Generalised_Arrhenius(A,n,E_a)
{
    // 1. Primary species data
    s_p_ = s_p;
    if ( s_p_ < 0.0 || s_p_ > 1.0 ) {
	ostringstream ost;
	ost << "Park_nonequilibrium::Park_nonequilibrium()" << endl
	    << "s_p was given as: " << s_p_ << " but is required to be between 0 and 1" << endl;
	input_error( ost );
    }
    iTp_ = iTp;
    
    // 2. Secondary species data
    iTq_ = iTq;
}

Park_nonequilibrium::
~Park_nonequilibrium() {}

int
Park_nonequilibrium::
s_eval(const Gas_data &Q)
{
    double T = pow( Q.T[iTp_], s_p_ ) * pow( Q.T[iTq_], 1.0 - s_p_ );
    
    Generalised_Arrhenius::eval_from_T(T);
    
    return SUCCESS;
}

Reaction_rate_coefficient* create_Park_nonequilibrium_coefficient(lua_State *L, Gas_model &g)
{
    return new Park_nonequilibrium(L, g);
}
