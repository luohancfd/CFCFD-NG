// Author: Daniel F. Potter
// Date: 10-December-2009

#include <cmath>
#include <sstream>
#include <iostream>


#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "gas-model.hh"
#include "GuptaYos-dcm.hh"
#include "physical_constants.hh"
#include "chemical-species-library.hh"

using namespace std;

GuptaYos_dcm::
GuptaYos_dcm(lua_State *L)
{
    // 1. Initialise species data from chemical species library
    if ( ! chemical_species_library_initialised() ) {
    	ostringstream ost;
	ost << "GuptaYos_dcm::GuptaYos_dcm()\n";
	ost << "The chemical species library must be initialised.\n";
	input_error(ost);
    }
    nsp_ = get_library_nsp();
    
    Z_.resize(nsp_);
    m_.resize(nsp_);
    x_.resize(nsp_);
    BI_table_.resize(nsp_);
    for ( int i = 0; i<nsp_; i++ ) BI_table_[i].resize(nsp_);

    for ( int isp = 0; isp < nsp_; ++isp ) {
    	Chemical_species * X = get_library_species_pointer(isp);

	// Get the species charge
	Z_[isp] = X->get_Z();
	
	// Get the particle mass
	m_[isp] = X->get_M() / PC_Avogadro;
    }
    

    lua_getglobal(L, "ignore_mole_fraction");
    if ( lua_isnil(L, -1) ) {
	cout << "GuptaYos_dcm::GuptaYos_dcm()\n";
	cout << "Setting ignore_mole_fraction_ to default value: "
	     << DEFAULT_MIN_MOLE_FRACTION << endl;
	ignore_mole_fraction_ = DEFAULT_MIN_MOLE_FRACTION;
    }
    else if ( lua_isnumber(L, -1) ) {
	ignore_mole_fraction_ = lua_tonumber(L, -1);
	if ( ignore_mole_fraction_ <= 0.0 ) {
	    ostringstream ost;
	    ost << "The ignore_mole_fraction is set to a value less than or equal to zero: " 
	        << ignore_mole_fraction_ << endl;
	    ost << "This must be a value between 0 and 1.\n";
	    input_error(ost); 
	}
	if ( ignore_mole_fraction_ > 1.0 ) {
	    ostringstream ost;
	    ost << "The ignore_mole_fraction is set to a value greater than 1.0: " 
	        << ignore_mole_fraction_ << endl;
	    ost << "This must be a value between 0 and 1.\n";
	    input_error(ost);
	}
    }
    else {
	ostringstream ost;
	ost << "The type ignore_mole_fraction is not a numeric type.\n";
	ost << "This must be a value between 0 and 1.\n";
	input_error(ost);
    }
    lua_pop(L, 1); // Pops "ignore_mole_fraction" off stack

    lua_getglobal(L, "collision_integrals");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "GuptaYos_dcm::GuptaYos_dcm()\n";
	ost << "Expected a table entry for 'collision_integrals'" << endl;
	input_error(ost);
    }
    
    for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
    	lua_rawgeti(L, -1, i);
    	string i_name = get_string(L, -1, "i");
    	string j_name = get_string(L, -1, "j");
    	int isp = get_library_index_from_name( i_name );
    	int jsp = get_library_index_from_name( j_name );
        unique_BIs_.push_back( new Binary_interaction(isp,jsp,L) );
	BI_table_[isp][jsp] = unique_BIs_.back();
	BI_table_[jsp][isp] = unique_BIs_.back();
	lua_pop(L, 1);	// pop collision_integral table item 'i'
    }
    lua_pop(L, 1);	// pop collision_integrals section
}

GuptaYos_dcm::
~GuptaYos_dcm()
{
    for ( size_t ic = 0; ic < unique_BIs_.size(); ++ic )
	delete unique_BIs_[ic];
}

int
GuptaYos_dcm::
s_eval_diffusion_coefficients(Gas_data &Q)
{
    // 0. Calculate mol-fractions
    convert_massf2molef( Q.massf, m_, x_ );
    
    // 1. Store uncorrected diffusion coefficients for all unique collision pairs
    for ( size_t ic=0; ic<unique_BIs_.size(); ++ic ) {
    	// Set D to zero if insufficient particles present
    	if ( x_[unique_BIs_[ic]->isp_] < ignore_mole_fraction_ || 
    	     x_[unique_BIs_[ic]->jsp_] < ignore_mole_fraction_ ) unique_BIs_[ic]->set_D(0.0);
    	else unique_BIs_[ic]->store_D(Q);
    }

    // 2. Map results onto the gas-data D_AB matrix
    for(int isp = 0; isp< nsp_; ++isp) {
	for(int jsp = 0; jsp< nsp_; ++jsp) {
	     Q.D_AB[isp][jsp] = BI_table_[isp][jsp]->get_D();
	}
    }

    return SUCCESS;
}

