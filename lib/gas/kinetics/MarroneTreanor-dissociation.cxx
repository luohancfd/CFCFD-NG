// Author: Daniel F Potter
// Date: 07-Dec-2009

#include <cmath>
#include <sstream>
#include <iostream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "MarroneTreanor-dissociation.hh"
#include "../models/physical_constants.hh"
#include "../models/chemical-species-library.hh"

using namespace std;

MarroneTreanor_dissociation::
MarroneTreanor_dissociation(lua_State *L, Gas_model &g)
    : Generalised_Arrhenius(L,g)
{
    U_ = get_number(L,-1,"U");

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
    	cout << "MarroneTreanor_dissociation::MarroneTreanor_dissociation()" << endl
    	     << "mode: vibration not found for species: "
    	     << X->get_name() << endl
    	     << "Exiting." << endl;
    	exit( BAD_INPUT_ERROR );
    }

    iTv_ = vib_modes_[0]->get_iT();
}

MarroneTreanor_dissociation::
MarroneTreanor_dissociation(double A, double n, double E_a, double U, string v_name)
    : Generalised_Arrhenius(A,n,E_a)
{
    U_ = U;
    
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
    	cout << "MarroneTreanor_dissociation::MarroneTreanor_dissociation()" << endl
    	     << "mode: vibration not found for species: "
    	     << X->get_name() << endl
    	     << "Exiting." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    iTv_ = vib_modes_[0]->get_iT();
}

MarroneTreanor_dissociation::
~MarroneTreanor_dissociation() {}

int
MarroneTreanor_dissociation::
s_eval(const Gas_data &Q)
{
    
    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    double T = Q.T[0];
    double Tv = Q.T[iTv_];
    
    // 1. Calculate pseudo-temperature gamma
    double gamma_inv = 1.0/Tv - 1.0/T - 1.0/U_;
    double gamma = 1.0 / gamma_inv;
    
    // 2. Calculate nonequilibrium factor
    double Z = 1.0;
    for ( size_t i=0; i<vib_modes_.size(); ++i ) {
    	Z *= vib_modes_[i]->eval_Q_from_T(T) * vib_modes_[i]->eval_Q_from_T(gamma) / 
    		vib_modes_[i]->eval_Q_from_T(Tv) / vib_modes_[i]->eval_Q_from_T(U_);
    }
    
    // 2. Evaluate GA coefficient
    Generalised_Arrhenius::eval_from_T(T);
    
    // 3. Augment k with NEQ factor
    k_ *= Z;
    
    return SUCCESS;
}

Reaction_rate_coefficient* create_MarroneTreanor_dissociation_coefficient(lua_State *L, Gas_model &g)
{
    return new MarroneTreanor_dissociation(L, g);
}
