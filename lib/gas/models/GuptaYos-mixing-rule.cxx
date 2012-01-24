// Author: Daniel F Potter
// Date: 13-Oct-2009

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "gas-model.hh"
#include "GuptaYos-mixing-rule.hh"
#include "physical_constants.hh"
#include "chemical-species-library.hh"

using namespace std;

GuptaYos_mixing_rule::
GuptaYos_mixing_rule( lua_State *L)
{
    // 1. Initialise species data from chemical species library
    if ( ! chemical_species_library_initialised() ) {
    	ostringstream ost;
	ost << "GuptaYos_mixing_rule::GuptaYos_mixing_rule()\n";
	ost << "The chemical species library must be initialised.\n";
	input_error(ost);
    }
    nsp_ = get_library_nsp();
    
    species_.resize(nsp_);
    m_.resize(nsp_, 0.0);
    Z_.resize(nsp_, 0 );
    x_.resize(nsp_, 0.0);
    BI_table_.resize(nsp_);
    for ( int i = 0; i<nsp_; i++ ) BI_table_[i].resize(nsp_);

    e_index_ = -1;
    for ( int isp = 0; isp < nsp_; ++isp ) {
    	Chemical_species * X = get_library_species_pointer(isp);
    	species_[isp] = get_library_species_pointer(isp);
	if ( X->get_name()=="e_minus" ) e_index_ = isp;

	// Get molecular weight and charge
	m_[isp] = X->get_M() / PC_Avogadro;	// convert kg/mol -> kg/particle
	Z_[isp] = X->get_Z();

	// Check that the species viscosity model is set to "collision integrals"
	lua_getglobal(L, X->get_name().c_str() );
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "GuptaYos_mixing_rule::GuptaYos_mixing_rule():\n";
	    ost << "Error locating information table for species: " << X->get_name() << endl;
	    input_error(ost);
	}
	lua_getfield(L, -1, "viscosity");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "GuptaYos_mixing_rule::GuptaYos_mixing_rule()\n";
	    ost << "Error setting viscosity table for species: " << X->get_name() << endl;
	    input_error(ost);
	}

	string vmodel = get_string(L, -1, "model");

	//cout << "vmodel= " << vmodel << endl;
	
	if ( vmodel != "collision integrals" ) {
	    ostringstream ost;
	    ost << "The viscosity model set for species: " << X->get_name() << endl;
	    ost << vmodel << " is not compatible with the Gupta-Yos mixing rule.\n";
	    input_error(ost);
	}
	lua_pop(L, 1); // Pop 'viscosity' table off stack

	// Check that the species conductivity model is set to "collision integrals"
	lua_getfield(L, -1, "thermal_conductivity");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "GuptaYos_mixing_rule::GuptaYos_mixing_rule()\n";
	    ost << "Error setting thermal_conductivity table for species: " << X->get_name() << endl;
	    input_error(ost);
	}

	string kmodel = get_string(L, -1, "model");
	
	if ( kmodel != "collision integrals" ) {
	    ostringstream ost;
	    ost << "The thermal conductivity model set for species: " << X->get_name() << endl;
	    ost << vmodel << " is not compatible with the Gupta-Yos mixing rule.\n";
	    input_error(ost);
	}
	lua_pop(L, 1); // Pop 'thermal conductivity' table off stack
	
	lua_pop(L, 1); // Pops "X" off stack
    }
    
    if ( e_index_ != -1 ) {
    	nsp_se_ = nsp_ - 1;
    	if ( e_index_ != nsp_ -1 ) {
    	    ostringstream ost;
    	    ost << "GuptaYos_mixing_rule::GuptaYos_mixing_rule()\n";
    	    ost << "Species 'e_minus' must be the last species in the list!" << endl;
    	    input_error(ost);
    	}
    }
    else {
    	nsp_se_ = nsp_;
    }
    lua_getglobal(L, "ignore_mole_fraction");
    if ( lua_isnil(L, -1) ) {
	cout << "GuptaYos_mixing_rule::GuptaYos_mixing_rule()\n";
	cout << "Setting ignore_mole_fraction to default value: "
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
	ost << "GuptaYos_mixing_rule::GuptaYos_mixing_rule()\n";
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

GuptaYos_mixing_rule::
~GuptaYos_mixing_rule()
{
    for ( size_t ic = 0; ic < unique_BIs_.size(); ++ic )
	delete unique_BIs_[ic];
}

Binary_interaction *
GuptaYos_mixing_rule::
get_binary_interaction_ptr( int isp, int jsp )
{
    for ( size_t i=0; i<unique_BIs_.size(); ++i ) {
    	if ( ( isp==unique_BIs_[i]->isp_ && jsp==unique_BIs_[i]->jsp_ ) ||
    	     ( isp==unique_BIs_[i]->jsp_ && jsp==unique_BIs_[i]->isp_ ) )
    	    return unique_BIs_[i];
    }
    
    cout << "GuptaYos_mixing_rule::get_binary_interaction_ptr()" << endl
         << "Species pair with indices: " << isp << " - " << jsp << " not found." << endl
         << "Exiting program." << endl;
    exit( BAD_INPUT_ERROR );
}

Collision_integral *
GuptaYos_mixing_rule::
get_collision_integral_ptr( int isp, int jsp )
{
    return get_binary_interaction_ptr(isp,jsp)->get_CI_model_ptr();
}

// alpha parameter for GuptaYos transport property calculation
#define eval_alpha() \
    ( 1.0 + (1.0-(m_[isp]/m_[jsp]))*(0.45-2.54*(m_[isp]/m_[jsp]))/pow(1.0+(m_[isp]/m_[jsp]),2) )

int
GuptaYos_mixing_rule::
s_eval_transport_coefficients(Gas_data &Q)
{    
    // 0. Set all transport parameters to zero
    Q.mu = 0.0;
    for ( size_t itm=0; itm<Q.k.size(); ++itm )	Q.k[itm] = 0.0;
    
    // 1. Calculate mol-fractions
    convert_massf2molef( Q.massf, m_, x_ );
    
    // 2. Store Delta_1 and Delta_2 values for all collision pairs
    for ( size_t ic=0; ic<unique_BIs_.size(); ++ic ) {
    	// Skip if insufficient particles present
    	if ( x_[unique_BIs_[ic]->isp_] < ignore_mole_fraction_ || 
    	     x_[unique_BIs_[ic]->jsp_] < ignore_mole_fraction_ ) continue;
    	unique_BIs_[ic]->store_Delta_1(Q);
    	unique_BIs_[ic]->store_Delta_2(Q);
    }
    
    double numerator=0.0, denominator=0.0;
    
    // 3. Calculate viscosity
    for ( int isp=0; isp<nsp_; ++isp ) {
    	if ( x_[isp] < ignore_mole_fraction_ ) continue;
    	numerator = m_[isp] * x_[isp];
    	denominator = 0.0;
    	for ( int jsp=0; jsp<nsp_; ++jsp ) {
    	    if ( x_[jsp] < ignore_mole_fraction_ ) continue;
    	    denominator += x_[jsp] * BI_table_[isp][jsp]->get_Delta_2();
    	}
    	Q.mu += numerator / denominator;
    }
    
    // 4. Calculate translational conductivities
    for ( int isp=0; isp<nsp_; ++isp ) {
    	if ( x_[isp] < ignore_mole_fraction_ ) continue;
    	numerator = x_[isp];
    	denominator = 0.0;
    	for ( int jsp=0; jsp<nsp_; ++jsp ) {
    	    if ( x_[jsp] < ignore_mole_fraction_ ) continue;
    	    denominator += eval_alpha() * x_[jsp] * BI_table_[isp][jsp]->get_Delta_2();
    	}
        // 15/4 = 3.75
    	Q.k[species_[isp]->get_iT_trans()] += ( numerator / denominator ) * PC_k_SI * 3.75;
    }
    
    // 5. Calculate (partially-excited) electronic conductivities
    for ( int isp=0; isp<nsp_; ++isp ) {
    	if ( x_[isp] < ignore_mole_fraction_ ) continue;
    	numerator = x_[isp] * species_[isp]->eval_Cv_elec(Q) / species_[isp]->get_R();
    	denominator = 0.0;
    	for ( int jsp=0; jsp<nsp_; ++jsp ) {
    	    if ( x_[jsp] < ignore_mole_fraction_ ) continue;
    	    denominator += x_[jsp] * BI_table_[isp][jsp]->get_Delta_1();
    	}
    	Q.k[species_[isp]->get_iT_elec()] += ( numerator / denominator ) * PC_k_SI;
    }
    
    // 6.  Calculate vibrational conductivities
    // 6a. Diatomics
    for ( int id=0; id<get_library_ndiatoms(); ++id ) {
    	Diatomic_species * XX = get_library_diatom_pointer(id);
    	int isp = XX->get_isp();
    	if ( x_[isp] < ignore_mole_fraction_ ) continue;
    	numerator = x_[isp] * XX->eval_Cv_vib(Q) / XX->get_R();
    	denominator = 0.0;
    	for ( int jsp=0; jsp<nsp_; ++jsp ) {
    	    if ( x_[jsp] < ignore_mole_fraction_ ) continue;
    	    denominator += x_[jsp] * BI_table_[isp][jsp]->get_Delta_1();
    	}
        Q.k[XX->get_iT_vib()] += ( numerator / denominator ) * PC_k_SI;
    }
    // 6b. Polyatomics
    for ( int ip=0; ip<get_library_npolyatoms(); ++ip ) {
    	Polyatomic_species * XXX = get_library_polyatom_pointer(ip);
    	int isp = XXX->get_isp();
    	if ( x_[isp] < ignore_mole_fraction_ ) continue;
    	numerator = x_[isp] * XXX->eval_Cv_vib(Q) / XXX->get_R();
    	denominator = 0.0;
    	for ( int jsp=0; jsp<nsp_; ++jsp ) {
    	    if ( x_[jsp] < ignore_mole_fraction_ ) continue;
    	    denominator += x_[jsp] * BI_table_[isp][jsp]->get_Delta_1();
    	}
        Q.k[XXX->get_iT_vib()] += ( numerator / denominator ) * PC_k_SI;
    }
    
    // 7. Calculate (fully-excited) rotational conductivities
    // 7a. Diatomics
    for ( int id=0; id<get_library_ndiatoms(); ++id ) {
    	Diatomic_species * XX = get_library_diatom_pointer(id);
    	int isp = XX->get_isp();
    	if ( x_[isp] < ignore_mole_fraction_ ) continue;
    	numerator = x_[isp];
    	denominator = 0.0;
    	for ( int jsp=0; jsp<nsp_; ++jsp ) {
    	    if ( x_[jsp] < ignore_mole_fraction_ ) continue;
    	    denominator += x_[jsp] * BI_table_[isp][jsp]->get_Delta_1();
    	}
    	Q.k[XX->get_iT_rot()] += ( numerator / denominator ) * PC_k_SI;
    }
    // 7b. Polyatomics
    for ( int ip=0; ip<get_library_npolyatoms(); ++ip ) {
    	Polyatomic_species * XXX = get_library_polyatom_pointer(ip);
    	int isp = XXX->get_isp();
    	if ( x_[isp] < ignore_mole_fraction_ ) continue;
    	numerator = x_[isp];
    	denominator = 0.0;
    	for ( int jsp=0; jsp<nsp_; ++jsp ) {
    	    if ( x_[jsp] < ignore_mole_fraction_ ) continue;
    	    denominator += x_[jsp] * BI_table_[isp][jsp]->get_Delta_1();
    	}
    	Q.k[XXX->get_iT_rot()] += ( numerator / denominator ) * PC_k_SI;
    }
    
    return SUCCESS;
}

GuptaYos_mixing_rule *
create_GuptaYos_mixing_rule_from_file( string lua_file )
{
    lua_State *L = initialise_lua_State();

    if( luaL_dofile(L, lua_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "GuptaYos_mixing_rule::create_GuptaYos_mixing_rule_from_file()\n";
	ost << "Error in gas model input file: " << lua_file << endl;
	input_error(ost);
    }
    
    return new GuptaYos_mixing_rule(L);
}
