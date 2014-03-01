// Author: Daniel F Potter
// Date: 07-Dec-2009
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "Knab-molecular-reaction.hh"
#include "../models/physical_constants.hh"
#include "../models/chemical-species-library.hh"

using namespace std;

Knab_molecular_reaction::
Knab_molecular_reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower)
    : Generalised_Arrhenius(L, g, T_upper, T_lower)
{
    U0_ = get_number(L, -1, "U0");
    U1_ = get_number(L, -1, "U1");
    alpha_ = get_number(L, -1, "alpha");
    alpha_A_ = alpha_*Generalised_Arrhenius::get_E_a()/PC_k_SI; // convert to K
    
    lua_getfield(L, -1, "Z_limit");
    if ( lua_isnumber(L, -1) ) {
	// User did provide Z_limit
	Z_limit_ = lua_tonumber(L, -1);
    }
    else {
	// User did not set Z_limit.
	// So set to -1.0 so it's not used.
	Z_limit_ = -1.0;
    }
    lua_pop(L, 1);

    string v_name = get_string(L,-1,"v_name");
    Chemical_species * X = get_library_species_pointer_from_name( v_name );
    // Search for the corresponding energy modes
    bool found = false;
    for ( int itm=0; itm<X->get_n_modes(); ++itm ) {
    	Species_energy_mode * sem = X->get_mode_pointer(itm);
    	if ( sem->get_type()=="vibration" ) {
    	    found = true;
    	    vib_modes_.push_back( sem );
    	}
    }
    if ( !found ) {
    	cout << "Knab_molecular_reaction::Knab_molecular_reaction()" << endl
    	     << "mode: vibration not found for species: "
    	     << X->get_name() << endl
    	     << "Exiting." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    iTv_ = vib_modes_[0]->get_iT();

    // 3. Set the reaction rate coefficient type
    type_ = "dissociation";
}

Knab_molecular_reaction::
Knab_molecular_reaction(double A, double n, double E_a, double T_upper, double T_lower,
			double U0, double U1, double alpha, double Z_limit, string v_name)
    : Generalised_Arrhenius(A, n, E_a, T_upper, T_lower), U0_(U0), U1_(U1), alpha_(alpha), alpha_A_(alpha*E_a/PC_k_SI), Z_limit_(Z_limit)
{
    Chemical_species * X = get_library_species_pointer_from_name( v_name );
    // Search for the corresponding energy modes
    bool found = false;
    for ( int itm=0; itm<X->get_n_modes(); ++itm ) {
    	Species_energy_mode * sem = X->get_mode_pointer(itm);
    	if ( sem->get_type()=="vibration" ) {
    	    found = true;
    	    vib_modes_.push_back( sem );
    	}
    }
    if ( !found ) {
    	cout << "Knab_molecular_reaction::Knab_molecular_reaction()" << endl
    	     << "mode: vibration not found for species: "
    	     << X->get_name() << endl
    	     << "Exiting." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    iTv_ = vib_modes_[0]->get_iT();

    // 3. Set the reaction rate coefficient type
    type_ = "dissociation";
}

Knab_molecular_reaction::
~Knab_molecular_reaction() {}

int
Knab_molecular_reaction::
s_eval(const Gas_data &Q)
{
    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    double T = Q.T[0];
    double Tv = Q.T[iTv_];
    
    // 1. Calculate pseudo-temperatures
    double U = U0_ + U1_*T;
    double gamma, T0, T_ast;
    if ( U > 0.0 ) {
	gamma = 1.0 / (1.0/Tv - 1.0/T - 1.0/U);
	T0 = 1.0 / (1.0/Tv - 1.0/U);
	T_ast = 1.0 / (1.0/T - 1.0/U);
    }
    else {
	// This is the special case of U=inf
	gamma = 1.0 / (1.0/Tv - 1.0/T);
	T0 = Tv;
	T_ast = T;
    }
    
    // 2. Calculate partition functions
    double Q_d_T = 1.0, Q_d_Tv = 1.0, Q_d_T0 = 1.0, Q_d_T_ast = 1.0;
    double Q_a_gamma = 1.0, Q_a_U = 1.0, Q_a_T0 = 1.0, Q_a_T_ast = 1.0;
    
    for ( size_t i=0; i<vib_modes_.size(); ++i ) {
    	Q_d_T *= vib_modes_[i]->eval_Q_from_T(T);
    	Q_d_Tv *= vib_modes_[i]->eval_Q_from_T(Tv);
    	Q_d_T0 *= vib_modes_[i]->eval_Q_from_T(T0);
    	Q_d_T_ast *= vib_modes_[i]->eval_Q_from_T(T_ast);
    	Q_a_gamma *= vib_modes_[i]->eval_Q_from_T(gamma,alpha_A_);
    	Q_a_T0 *= vib_modes_[i]->eval_Q_from_T(T0,alpha_A_);
    	Q_a_T_ast *= vib_modes_[i]->eval_Q_from_T(T_ast,alpha_A_);
	if ( U > 0.0 ) {
	    Q_a_U *= vib_modes_[i]->eval_Q_from_T(-U,alpha_A_);
	}
	else {
	    // Special case of U = inf
	    // I showed this numerically by plotting
	    // partition function of truncated harmonic
	    // oscillator for large values of T
	    Q_a_U *= alpha_A_/vib_modes_[i]->get_theta();
	}
    }
    
    // 3. Calculate nonequilibrium factor
    double tmpA = Q_d_T / Q_d_Tv;
    double tmpB = exp(-alpha_A_/T)*Q_a_gamma + Q_d_T0 - Q_a_T0;
    double tmpC = exp(-alpha_A_/T)*Q_a_U + Q_d_T_ast - Q_a_T_ast;
    double Z = tmpA * tmpB / tmpC;
    // At small values of temperature, the values for
    // Q_d_T0 and Q_a_T0 are essentially the same (to machine precision).
    // Similarly, are Q_d_T_ast and Q_a_T_ast essentially the same. Also, the
    // exponential term is very small (<1.0e-50).  All of these factors conspire
    // to give numer = 0.0 and denom = 0.0 which leads to Z = nan.
    // The fix is to set Z = 1.0 in these instances.  This will be reasonable
    // because this occure at low temperatures when T does not differ greatly from Tvib.
    // In these case the deviation form equiilibrium won't  be very strong and so
    // a value of Z = 1.0 will be ok.

    if( std::isnan(Z) || std::isinf(Z) )
	Z = 1.0;

    // 4. Evaluate GA coefficient
    // Check on temperature limits
    if ( T > T_upper_ )
	T = T_upper_;
    if ( T < T_lower_ )
	T = T_lower_;
    Generalised_Arrhenius::eval_from_T(T);
    // 5. Augment k with NEQ factor
    if ( Z_limit_ > 0.0 && Z > Z_limit_ ) {
	// We put a limit on the nonequilibrium factor in case it
	// gets too large (that is, crazy) for some flow conditions
	Z = Z_limit_;
    }
    k_ *= fabs(Z);

    return SUCCESS;
}

Reaction_rate_coefficient* create_Knab_molecular_reaction_coefficient(lua_State *L, Gas_model &g, double T_upper, double T_lower)
{
    return new Knab_molecular_reaction(L, g, T_upper, T_lower);
}
