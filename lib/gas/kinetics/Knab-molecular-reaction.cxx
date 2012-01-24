// Author: Daniel F Potter
// Date: 07-Dec-2009

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
Knab_molecular_reaction(lua_State *L, Gas_model &g)
    : Generalised_Arrhenius(L,g)
{
}

Knab_molecular_reaction::
Knab_molecular_reaction(double A, double n, double E_a, double U, double alpha, double A_var, string v_name)
    : Generalised_Arrhenius(A,n,E_a)
{
    U_ = U;
    alpha_ = alpha;
    A_var_ = A_var;
    
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
    double gamma = 1.0 / ( 1.0/Tv - 1.0/T - 1.0/U_ );
    double T0 = 1.0 / ( 1.0/Tv - 1.0/U_ );
    double T_ast = 1.0 / ( 1.0/T - 1.0/U_ );
    
    // 2. Calculate partition functions
    double Q_d_T = 1.0, Q_d_Tv = 1.0, Q_d_T0 = 1.0, Q_d_T_ast = 1.0;
    double Q_a_gamma = 1.0, Q_a_U = 1.0, Q_a_T0 = 1.0, Q_a_T_ast = 1.0;
    
    for ( size_t i=0; i<vib_modes_.size(); ++i ) {
    	Q_d_T *= vib_modes_[i]->eval_Q_from_T(T);
    	Q_d_Tv *= vib_modes_[i]->eval_Q_from_T(Tv);
    	Q_d_T0 *= vib_modes_[i]->eval_Q_from_T(T0);
    	Q_d_T_ast *= vib_modes_[i]->eval_Q_from_T(T_ast);
    	Q_a_gamma *= vib_modes_[i]->eval_Q_from_T(gamma);
    	Q_a_U *= vib_modes_[i]->eval_Q_from_T(-U_,alpha_*A_var_);
    	Q_a_T0 *= vib_modes_[i]->eval_Q_from_T(T0,alpha_*A_var_);
    	Q_a_T_ast *= vib_modes_[i]->eval_Q_from_T(T_ast,alpha_*A_var_);
    }
    
    // 3. Calculate nonequilibrium factor
    double tmpA = Q_d_T / Q_d_Tv;
    double tmpB = exp( - alpha_ * A_var_ / T ) * Q_a_gamma + Q_d_T0 - Q_a_T0;
    double tmpC = exp( - alpha_ * A_var_ / T ) * Q_a_U + Q_d_T_ast - Q_a_T_ast;
    double Z = tmpA * tmpB / tmpC;
    
    // 4. Evaluate GA coefficient
    Generalised_Arrhenius::eval_from_T(T);
    
    // 3. Augment k with NEQ factor
    k_ *= fabs(Z);
    
    return SUCCESS;
}

Reaction_rate_coefficient* create_Knab_molecular_reaction_coefficient(lua_State *L, Gas_model &g)
{
    return new Knab_molecular_reaction(L, g);
}
