// Author: Rowan J. Gollan
// Date: 12-Sep-2008

#include <iostream>
#include <sstream>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"
#include "../../nm/source/no_fuss_linear_algebra.hh"

#include "../models/chemical-species-library.hh"
#include "../models/chemical-species.hh"

#include "chemical-kinetic-MC-system.hh"
#include "chemistry-energy-coupling.hh"

#define CHECK_GAS_DATA 0
#define APPLY_AVERAGE_ENERGY_SOURCE_TERMS 1

using namespace std;

Chemical_kinetic_MC_system::
Chemical_kinetic_MC_system(lua_State *L, Gas_model &g, int nreac, double error_tol)
    : OdeSystem(2*nreac, true), err_tol_(error_tol)
{
    lua_getglobal(L, "reactions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Chemical_kinetic_MC_system::Chemical_kinetic_MC_system()\n";
	ost << "Error interpreting 'reactions'; a table of reactions is expected.\n";
	input_error(ost);
    }

    vector<Coupling_component*> ccs;
    
    for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i);
	reaction_.push_back(create_Reaction(L, g));
	// NOTE: 'i-1' is required due to lua +1 offset
	create_Coupling_components_for_reaction( L, reaction_.back(), int(i-1), ccs );
	lua_pop(L, 1);
    }
    
    cout << "Chemical_kinetic_MC_system::Chemical_kinetic_MC_system()" << endl
         << "-> found " << ccs.size() << " Coupling_component's" << endl;
    
    lua_pop(L, 1);	// pop reactions

    int nsp = g.get_number_of_species();
    participation_.resize(nsp);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
    	vector<int> nu;
	for ( int ir = 0; ir < nreac; ++ir ) {
	    int nu_test = reaction_[ir]->get_nu(isp);
	    if ( nu_test != 0 ) {
		participation_[isp].push_back(ir);
		nu.push_back( nu_test );
	    }
	}
	// Now we can create the SpeciesPieces
	spec_.push_back( new Species_pieces( isp, participation_[isp], nu ) );
    }
    
    // Group the Coupling_components into Chemistry_energy_coupling mechanisms,
    // one for each thermal mode of each species
    if ( ccs.size()>0 ) {
    	for ( int isp=0; isp<nsp; ++isp ) {
    	    Chemical_species * X = get_library_species_pointer(isp);
    	    for ( int itm=0; itm<X->get_n_modes(); ++itm ) {
    	    	string mode = X->get_mode_pointer(itm)->get_type();
    	    	create_Chemistry_energy_coupling_for_species_mode(isp,mode,ccs,cecs_);
    	    }
    	}
    }
    
    cout << "Chemical_kinetic_MC_system::Chemical_kinetic_MC_system()" << endl
         << "-> created " << cecs_.size() << " Chemistry_energy_coupling's" << endl;

    ydot_.resize(2*nreac, 0.0);
    w_.resize(nreac, 0.0);
    cinit_.resize(nsp, 0.0);
    c_.resize(nsp, 0.0);
    massf_.resize(nsp, 0.0);
    M_.resize(nsp, 0.0);
    for ( size_t isp = 0; isp < M_.size(); ++isp ) M_[isp] = g.molecular_weight(isp);
}

Chemical_kinetic_MC_system::
Chemical_kinetic_MC_system(string cfile, Gas_model &g, int nreac, double error_tol)
    : OdeSystem(2*nreac, true), err_tol_(error_tol)
{
    // Do some pre-work on the cfile to massage
    // it into a state to be parsed for the
    // individual reactions.

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    // Set up species table
    lua_newtable(L);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	// This species table maps to C++ indices, because
	// it is used to setup the integer maps for
	// the reaction coefficients.
	lua_pushinteger(L, isp);
	lua_setfield(L, -2, g.species_name(isp).c_str());
    }
    // Plus add a field 'size': no of species
    lua_pushinteger(L, g.get_number_of_species());
    lua_setfield(L, -2, "size");
    lua_setglobal(L, "species");

    // Path to reaction parsing script
    string home(getenv("HOME"));
    string script_file(home);
    script_file.append("/e3bin/reaction_parser.lua");

    if ( luaL_dofile(L, script_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Chemical_kinetic_MC_system():\n";
	ost << "Error in loading script file: " << script_file << endl;
	input_error(ost);
    }
    
    // Parse the input file...
    lua_getglobal(L, "main");
    lua_pushstring(L, cfile.c_str());
    if ( lua_pcall(L, 1, 0, 0) != 0 ) {
	ostringstream ost;
	ost << "Chemical_kinetic_MC_system():\n";
	ost << "Error trying to load reaction scheme file: " << cfile << endl;
	ost << "Lua error message: " << lua_tostring(L, -1) << endl;
	input_error(ost);
    }
    

    lua_getglobal(L, "reactions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Chemical_kinetic_MC_system::Chemical_kinetic_MC_system()\n";
	ost << "Error interpreting 'reactions'; a table of reactions is expected.\n";
	input_error(ost);
    }
    
    vector<Coupling_component*> ccs;
    
    for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i);
	reaction_.push_back(create_Reaction(L, g));
	create_Coupling_components_for_reaction( L, reaction_.back(), int(i), ccs );
	lua_pop(L, 1);
    }
    lua_pop(L, 1);
    
    int nsp = g.get_number_of_species();
    participation_.resize(nsp);
    for ( int isp = 0; isp < nsp; ++isp ) {
    	vector<int> nu;
	for ( int ir = 0; ir < nreac; ++ir ) {
	    int nu_test = reaction_[ir]->get_nu(isp);
	    if ( nu_test != 0 ) {
		participation_[isp].push_back(ir);
		nu.push_back( nu_test );
	    }
	}
	// Now we can create the SpeciesPieces
	spec_.push_back( new Species_pieces( isp, participation_[isp], nu ) );
    }
    
    // Group the Coupling_components into Chemistry_energy_coupling mechanisms,
    // one for each thermal mode of each species
    if ( ccs.size()>0 ) {
    	for ( int isp=0; isp<nsp; ++isp ) {
    	    Chemical_species * X = get_library_species_pointer(isp);
    	    for ( int itm=0; itm<X->get_n_modes(); ++itm ) {
    	    	string mode = X->get_mode_pointer(itm)->get_type();
    	    	create_Chemistry_energy_coupling_for_species_mode(isp,mode,ccs,cecs_);
    	    }
    	}
    }
    
    cout << "Chemical_kinetic_MC_system::Chemical_kinetic_MC_system()" << endl
         << "-> created " << cecs_.size() << " Chemistry_energy_coupling's" << endl;

    ydot_.resize(2*nreac, 0.0);
    w_.resize(nreac, 0.0);
    cinit_.resize(nsp, 0.0);
    c_.resize(nsp, 0.0);
    massf_.resize(nsp, 0.0);
    M_.resize(nsp, 0.0);
    for ( size_t isp = 0; isp < M_.size(); ++isp ) M_[isp] = g.molecular_weight(isp);
    
    lua_close(L);
}

Chemical_kinetic_MC_system::
~Chemical_kinetic_MC_system()
{
    for ( size_t i = 0; i < reaction_.size(); ++i ) {
	delete reaction_[i];
    }
    
    for ( size_t i = 0; i < spec_.size(); ++i ) {
	delete spec_[i];
    }
    
    for ( size_t i = 0; i < cecs_.size(); ++i ) {
    	delete cecs_[i];
    }
}

int
Chemical_kinetic_MC_system::
eval(const valarray<double> &y, valarray<double> &ydot)
{
    // Firstly compute the rate coefficients
    if ( ! called_at_least_once ) {
	for ( size_t ir = 0; ir < reaction_.size(); ++ir ) {
	    if ( reaction_[ir]->compute_kf_first() ) {
		reaction_[ir]->compute_k_f(*Q_);
		reaction_[ir]->compute_k_b(*Q_);
	    }
	    else {
		reaction_[ir]->compute_k_b(*Q_);
		reaction_[ir]->compute_k_f(*Q_);
	    }
	}
	called_at_least_once = true;
    }
    
    // Now calculate concentration vector
    this->eval_new_concentrations( y, c_ );
    
    // Finally assemble the rate.
    for ( size_t ir = 0; ir < reaction_.size(); ++ir ) {
    	reaction_[ir]->compute_forward_rate( c_ );
    	reaction_[ir]->compute_backward_rate( c_ );
	ydot[2*ir] = reaction_[ir]->w_f();
	ydot[2*ir+1] = reaction_[ir]->w_b();
    }
    
    return SUCCESS;
}

int
Chemical_kinetic_MC_system::
eval_new_concentrations( const valarray<double> &y, valarray<double> &c )
{
    for ( size_t ir=0; ir<w_.size(); ++ir ) {
    	w_[ir] = y[2*ir] - y[2*ir+1];
    }
    
    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
	c[isp] = spec_[isp]->eval_conc( w_ );
    }
    
    return SUCCESS;
}

int
Chemical_kinetic_MC_system::
eval_new_concentrations( const valarray<double> &y, vector<double> &c )
{
    for ( size_t ir=0; ir<w_.size(); ++ir ) {
    	w_[ir] = y[2*ir] - y[2*ir+1];
    }
    
    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
	c[isp] = spec_[isp]->eval_conc( w_ );
    }
    
    return SUCCESS;
}

int
Chemical_kinetic_MC_system::
eval_species_rates( const valarray<double> &y, valarray<double> &cdot )
{
    eval(y,ydot_);
    int ir;
    
    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
    	cdot[isp] = 0.0;
	for( size_t i = 0; i < spec_[isp]->nu_.size(); ++i ) {
	    ir = spec_[isp]->reac_index_[i];
	    cdot[isp] += spec_[isp]->nu_[i] * ( ydot_[2*ir] - ydot_[2*ir + 1] );
	}
    }
    
    return SUCCESS;
}

const double eps1 = 0.001;
const double chem_step_upper_limit = 1.0e-3;
const double chem_step_lower_limit = 1.0e-20;
const double zero_tol = 1.0e-30;

double
Chemical_kinetic_MC_system::
stepsize_select(const valarray<double> &y)
{
    eval(y, ydot_);
    double min_dt = chem_step_upper_limit; // to get us started
    double old_dt = 0.0;
    double rate = 0.0;
    int ir;
    
    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
    	rate = 0.0;
	for( size_t i = 0; i < spec_[isp]->nu_.size(); ++i ) {
	    ir = spec_[isp]->reac_index_[i];
	    rate += spec_[isp]->nu_[i] * ( ydot_[2*ir] - ydot_[2*ir + 1] );
	}
    	
	if( (spec_[isp]->get_init_conc() > 0.0) && (fabs(rate) > zero_tol) ) {
	    old_dt = fabs( spec_[isp]->get_init_conc() / rate );
	    if( old_dt < min_dt ) {
		min_dt = old_dt;
	    }
	}
    }

    double dt_chem = eps1 * min_dt;

    // Impose upper and lower chem_step limits
    if( dt_chem > chem_step_upper_limit )
	dt_chem = chem_step_upper_limit;

    if ( dt_chem < chem_step_lower_limit )
	dt_chem = chem_step_lower_limit;

    return dt_chem;
}

const double min_conc = 1.0e-30;

bool
Chemical_kinetic_MC_system::
passes_system_test(valarray<double> &y)
{
    for ( size_t iy=0; iy<y.size(); ++iy ) {
    	if ( !finite( y[iy] ) ) return false;
    }
    return true;
}

void
Chemical_kinetic_MC_system::
print_reaction_rates()
{
    for ( size_t ir = 0; ir < reaction_.size(); ++ir ) {
	printf("rxn[%i]: w_f_=%16.15e, w_b_=%16.15e\n", int(ir), reaction_[ir]->w_f(), reaction_[ir]->w_b());
    }
    return;
}

void
Chemical_kinetic_MC_system::
print_species_rates()
{
    int ir;
    
    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
    	double rate = 0.0;
	for( size_t i = 0; i < spec_[isp]->nu_.size(); ++i ) {
	    ir = spec_[isp]->reac_index_[i];
	    rate += spec_[isp]->nu_[i] * ( ydot_[2*ir] - ydot_[2*ir + 1] );
	}
	double dt = fabs( spec_[isp]->get_init_conc() / rate );
	cout << "isp = " << isp << ", rate = " << rate << ", dt = " << dt << endl;
    }
    
    return;
}

void
Chemical_kinetic_MC_system::
print_limiting_species_and_reaction()
{
    double min_dt = chem_step_upper_limit;
    int isp_ast = -1, ir_ast = -1;
    int ir;
    
    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
	for( size_t i = 0; i < spec_[isp]->nu_.size(); ++i ) {
	    ir = spec_[isp]->reac_index_[i];
	    double rate = spec_[isp]->nu_[i] * ( ydot_[2*ir] - ydot_[2*ir + 1] );
	    if( (spec_[isp]->get_init_conc() > 0.0) && (fabs(rate) > zero_tol) ) {
	    	double dt = fabs( spec_[isp]->get_init_conc() / rate );
	    	if ( dt < min_dt ) {
	    	    min_dt = dt;
	    	    isp_ast = isp; ir_ast = ir;
	    	}
	    }
	}
    }
    
    cout << "limiting species and reaction: isp = " << isp_ast << ", ir_ast = " << ir_ast << ", min_dt = " << min_dt << endl;

    return;
}

void
Chemical_kinetic_MC_system::
set_gas_data_ptr_and_initial_concs( Gas_data &Q, vector<double> &c )
{
    // Set the gas data pointer
    Q_ = &Q;
    
    // Set the concentrations and map onto the species pieces
    for ( size_t isp=0; isp<spec_.size(); ++isp ) {
    	cinit_[isp] = c[isp];
    	spec_[isp]->set_init_conc( cinit_[isp] );
    }
    
    return;
}

int
Chemical_kinetic_MC_system::
initialise_chemistry_energy_coupling( Gas_data &Q, vector<double> &c_old )
{
    for ( size_t i=0; i<cecs_.size(); ++i ) {
    	cecs_[i]->set_e_and_N_old(Q,c_old);
    }
    
    return SUCCESS;
}

int
Chemical_kinetic_MC_system::
apply_chemistry_energy_coupling( Gas_data &Q, valarray<double> &delta_c, 
    				 vector<double> &c_new )
{
#   if CHECK_GAS_DATA
    if ( !Q.check_values() ) {
    	cout << "Gas data is bad prior to applying chemistry energy coupling." << endl;
    	return FAILURE;
    }
#   endif
    
    for ( size_t i=0; i<cecs_.size(); ++i ) {
    	if ( cecs_[i]->update_energy(Q,delta_c,c_new)==FAILURE )
    	    return FAILURE;
#       if CHECK_GAS_DATA
	if ( !Q.check_values() ) {
	    cout << "Gas data is bad after applying chemistry energy coupling " 
	         << "for species index: " << cecs_[i]->get_isp() << endl;
	    return FAILURE;
	}
#       endif
    }
    
    return SUCCESS;
}

int
Chemical_kinetic_MC_system::
eval_chemistry_energy_coupling_source_terms( Gas_data &Q, const valarray<double> &y, vector<double> &dedt )
{
    // NOTE: really don't need this eval call as eval_species_rates should have just been run
    eval(y,ydot_);
    
    int imode;
    for ( size_t i=0; i<cecs_.size(); ++i ) {
    	imode = cecs_[i]->get_imode();
    	dedt[imode] += cecs_[i]->eval_source_term(Q,ydot_);
    }
    
#   if APPLY_AVERAGE_ENERGY_SOURCE_TERMS
    // Now apply the source terms due to species taking the average modal 
    // energies with them as they are created or destroyed
    Chemical_species * X;
    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
	for( size_t i = 0; i < spec_[isp]->nu_.size(); ++i ) {
	    int ir = spec_[isp]->reac_index_[i];
	    double rate = spec_[isp]->nu_[i] * ( ydot_[2*ir] - ydot_[2*ir + 1] );
	    X = get_library_species_pointer(isp);
	    for ( int j=0; j<X->get_n_modes(); ++j ) {
	    	imode = X->get_mode_pointer(j)->get_iT();
	    	// Don't need to worry about the translational mode
	    	if ( imode>0 ) {
	    	    dedt[imode] += rate * X->get_mode_pointer(j)->eval_energy(Q) * X->get_M() / Q.rho;	// mole/m3-s * J/kg * kg/mole / kg/m3= W/kg
	    	    // cout << "dedt = " << rate * X->get_mode_pointer(j)->eval_energy(Q) * X->get_M() << endl;
	    	}
	    }
	}
    }
#   endif
    
    return SUCCESS;
}
