/** \file composite-gas-model.cxx
 *  \ingroup gas
 *
 *  \author Rowan J Gollan
 *  \version 09-Jul-2008
 **/

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "composite-gas-model.hh"
#include "perfect-gas-EOS.hh"
#include "simple-gas-EOS.hh"
#include "noble-abel-gas-EOS.hh"
#include "van-der-waals-gas-EOS.hh"
#include "Bender-EOS.hh"
#include "MBWR-EOS.hh"
#include "noneq-gas-EOS.hh"
#include "constant-specific-heats.hh"
#include "perfect-thermal-behaviour.hh"
#include "real-thermal-behaviour.hh"
#include "dense-real-thermal-behaviour.hh"
#include "noneq-thermal-behaviour.hh"
#include "Wilke-mixing-rule.hh"
#include "no-diffusion-coefficients.hh"
#include "hard-sphere-dcm.hh"
#include "sound-speed-model.hh"
#include "no-transport-coefficients.hh"
#include "GuptaYos-mixing-rule.hh"
#include "GuptaYos-dcm.hh"
#include "Armaly-Sutton-mixing-rule.hh"

using namespace std;

Composite_gas_model::
Composite_gas_model(string cfile)
    : Gas_model(cfile)
{
    lua_State *L = initialise_lua_State();
    
    if( luaL_dofile(L, cfile.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Composite_gas_model():\n";
	ost << "Error in gas model input file: " << cfile << endl;
	input_error(ost);
    }

    lua_getglobal(L, "species");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Composite_gas_model::Composite_gas_model():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }
    int nsp = lua_objlen(L, -1);
    lua_pop(L,1);	// pop 'species' off the stack
    set_number_of_species(nsp);

    string EOS = get_string(L, LUA_GLOBALSINDEX, "equation_of_state");
    if ( EOS == "perfect gas" ) {
	EOS_ = (Equation_of_state*) new Perfect_gas(L);
    }
    else if ( EOS == "Noble-Abel gas" ) {
	EOS_ = (Equation_of_state*) new Noble_Abel_gas(L);
    }
    else if ( EOS == "van der Waals gas" ) {
     	EOS_ = (Equation_of_state*) new van_der_Waals_gas(L);
    }
    else if ( EOS == "simple gas" ) {
        EOS_ = (Equation_of_state*) new Simple_gas(L);
    }
    else if ( EOS == "Bender" ) {
        EOS_ = (Equation_of_state*) new Bender_EOS(L);
    }
    else if ( EOS == "MBWR" ) {
        EOS_ = (Equation_of_state*) new MBWR_EOS(L);
    }
    else if ( EOS == "nonequilibrium gas" ) {
    	EOS_ = (Equation_of_state*) new Noneq_gas(L);
    }
    else {
	ostringstream ost;
	ost << "Composite_gas_model::Composite_gas_model():\n";
	ost << "The 'equation_of_state': " << EOS << endl;
	ost << "is unknown or not implemented.\n";
	input_error(ost);
    }
    
    string TBM = get_string(L, LUA_GLOBALSINDEX, "thermal_behaviour");
    if ( TBM == "constant specific heats" ) {
	TBM_ = (Thermal_behaviour_model*) new Constant_specific_heats(L);
    }
    else if ( TBM == "thermally perfect" ) {
	set_reaction_compatibility(true);
	TBM_ = (Thermal_behaviour_model*) new Perfect_thermal_behaviour(L);
    }
    else if ( TBM == "thermally real" ) {
	TBM_ = (Thermal_behaviour_model*) new Real_thermal_behaviour(L);
    }
    else if ( TBM == "dense thermally real" ) {
        TBM_ = (Thermal_behaviour_model*) new Dense_real_thermal_behaviour(L);
    }
    else if ( TBM == "thermal nonequilibrium" ) {
	set_reaction_compatibility(true);
	TBM_ = (Thermal_behaviour_model*) new Noneq_thermal_behaviour(L);
	Noneq_thermal_behaviour *ntb = dynamic_cast<Noneq_thermal_behaviour*>(TBM_);
	m_components_.resize(ntb->get_number_of_modes());
	for ( int imode = 0; imode < ntb->get_number_of_modes(); ++imode ) {
	    m_names_.push_back(ntb->mode_name(imode));
	    for ( int ic = 0; ic < ntb->mode_no_components(imode); ++ic ) {
		m_components_[imode].push_back(ntb->mode_component_name(imode, ic));
	    }
	}
    }
    else {
	ostringstream ost;
	ost << "Composite_gas_model::Composite_gas_model():\n";
	ost << "The 'thermal_behaviour': " << TBM << endl;
	ost << "is unknown or not implemented.\n";
	input_error(ost);
    }
    set_number_of_modes(TBM_->get_number_of_modes());

    string TCM = get_string(L, LUA_GLOBALSINDEX, "mixing_rule");
    if ( TCM == "Wilke" ) {
	TCM_ = (Transport_coefficients_model*) new Wilke_mixing_rule(L);
    }
    else if ( TCM == "GuptaYos" ) {
	TCM_ = (Transport_coefficients_model*) new GuptaYos_mixing_rule(L);
    }
    else if ( TCM == "ArmalySutton") {
	TCM_ = (Transport_coefficients_model*) new Armaly_Sutton_mixing_rule(L);
    }
    else if ( TCM == "None" ) {
    	TCM_ = (Transport_coefficients_model*) new No_transport_coefficients();
    }
    else {
	ostringstream ost;
	ost << "Composite_gas_model::Composite_gas_model():\n";
	ost << "The 'mixing_rule': " << TCM << endl;
	ost << "is unknown or not implemented.\n";
	input_error(ost);
    }

    if ( nsp > 1 ) {
	string DCM = get_string(L, LUA_GLOBALSINDEX, "diffusion_coefficients");
	if ( DCM == "hard sphere" ) {
	    DCM_ = (Diffusion_coefficients_model*) new Hard_sphere_dcm(L);
	}
	else if ( DCM == "GuptaYos" ) {
	    // check that transport model has also been set to GuptaYos
	    if ( TCM != "GuptaYos" ) {
		ostringstream ost;
		ost << "Composite_gas_model::Composite_gas_model():\n";
		ost << "For physical consistency, the transport model should also be set to 'GuptaYos'\n";
		ost << "for 'GuptaYos' diffusion coefficients.\n";
		input_error(ost);
	    }
	    DCM_ = (Diffusion_coefficients_model*) new GuptaYos_dcm(L);
	}
	else if ( DCM == "None" ) {
	    DCM_ = (Diffusion_coefficients_model*) new No_diffusion_coefficients();
	}
	else {
	    ostringstream ost;
	    ost << "Composite_gas_model::Composite_gas_model():\n";
	    ost << "The 'diffusion_coefficients' model: " << DCM << endl;
	    ost << "is unknown or not implemented.\n";
	    input_error(ost);
	}
    }
    else
	DCM_ = (Diffusion_coefficients_model*) new No_diffusion_coefficients();
    
    string SSM = get_string(L, LUA_GLOBALSINDEX, "sound_speed");
    SSM_ = create_sound_speed_model(SSM, *this, L);
    lua_close(L);
}

Composite_gas_model::
~Composite_gas_model()
{
    delete EOS_;
    delete TBM_;
    delete TCM_;
    delete DCM_;
    delete SSM_;
}

int
Composite_gas_model::
s_eval_thermo_state_rhoe(Gas_data &Q)
{
    // 0. Check we have useful inputs
    if ( Q.rho <= 0.0 )
	return FAILURE;

    // 1. Evalaute temperatures from energies
    if ( TBM_->eval_temperature(Q, EOS_) != SUCCESS )
	return FAILURE;

    // 2. With temperature known, call the equation of state
    if ( EOS_->eval_pressure(Q) != SUCCESS )
	return FAILURE;

    // 3. It is now possible to compute sound speed.
    if ( SSM_->eval_sound_speed(Q) != SUCCESS )
	return FAILURE;

    // if we got this far, assume everything was OK
    return SUCCESS;
}

int
Composite_gas_model::
s_eval_thermo_state_pT(Gas_data &Q)
{
    // 0. Check we have useful inputs
    if ( Q.p <= 0.0 )
	return FAILURE;

    if ( Q.T[0] <= 0.0 )
	return FAILURE;

    // 1. Evaluate the density from p and T
    if ( EOS_->eval_density(Q) != SUCCESS )
	return FAILURE;
    
    // 2. Compute the specific energies
    if ( TBM_->eval_energy(Q, EOS_) != SUCCESS )
	return FAILURE;
    
    // 3. Compute sound speed.
    if ( SSM_->eval_sound_speed(Q) != SUCCESS )
	return FAILURE;

    // if we got this far, assume everything was OK
    return SUCCESS;

}

int
Composite_gas_model::
s_eval_thermo_state_rhoT(Gas_data &Q)
{
    // 0. Check we have useful inputs
    if ( Q.rho <= 0.0 )
	return FAILURE;

    if ( Q.T[0] <= 0.0 )
	return FAILURE;

    // 1. Evaluate the pressure from rho and T
    if ( EOS_->eval_pressure(Q) != SUCCESS )
	return FAILURE;

    // 2. Compute the specific energies
    if ( TBM_->eval_energy(Q, EOS_) != SUCCESS )
	return FAILURE;

    // 3. Compute sound speed.
    if ( SSM_->eval_sound_speed(Q) != SUCCESS )
	return FAILURE;

    // if we got this far, assume everything was OK
    return SUCCESS;
}

int
Composite_gas_model::
s_eval_thermo_state_rhop(Gas_data &Q)
{
    // 0. Check we have useful inputs
    if ( Q.rho <= 0.0 )
	return FAILURE;

    if ( Q.p <= 0.0 )
	return FAILURE;

    // 1. Evaluate the temperature from rho and p
    if ( EOS_->eval_temperature(Q) != SUCCESS )
	return FAILURE;

    // 2. Compute the specific energies
    if ( TBM_->eval_energy(Q, EOS_) != SUCCESS )
	return FAILURE;

    // 3. Compute sound speed.
    if ( SSM_->eval_sound_speed(Q) != SUCCESS )
	return FAILURE;

    // if we got this far, assume everything was OK
    return SUCCESS;
}

int
Composite_gas_model::
s_eval_sound_speed(Gas_data &Q)
{
    return SSM_->eval_sound_speed(Q);
}

int
Composite_gas_model::
s_eval_transport_coefficients(Gas_data &Q, Gas_model *gmodel)
{
    return TCM_->eval_transport_coefficients(Q, gmodel);
}

int
Composite_gas_model::
s_eval_diffusion_coefficients(Gas_data &Q)
{
    return DCM_->eval_diffusion_coefficients(Q);
}

double
Composite_gas_model::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    return EOS_->dTdp_const_rho(Q, status);
}

double
Composite_gas_model::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    return EOS_->dTdrho_const_p(Q, status);
}

double
Composite_gas_model::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    return EOS_->dpdrho_const_T(Q, status);
}

double
Composite_gas_model::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    return EOS_->dpdrho_i_const_T(Q, isp, status);
}

double
Composite_gas_model::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    return EOS_->dpdT_i_const_rho(Q, itm, status);
}

double
Composite_gas_model::
s_dedT_const_v(const Gas_data &Q, int &status)
{
    return TBM_->dedT_const_v(Q, EOS_, status);
}

double
Composite_gas_model::
s_dhdT_const_p(const Gas_data &Q, int &status)
{
    return TBM_->dhdT_const_p(Q, EOS_, status);
}

double
Composite_gas_model::
s_internal_energy(const Gas_data &Q, int isp)
{
    return TBM_->eval_energy_isp(Q, EOS_, isp);
}

double
Composite_gas_model::
s_enthalpy(const Gas_data &Q, int isp)
{
    return TBM_->eval_enthalpy_isp(Q, EOS_, isp);
}

double
Composite_gas_model::
s_entropy(const Gas_data &Q, int isp)
{
    return TBM_->eval_entropy_isp(Q, EOS_, isp);
}

double
Composite_gas_model::
s_modal_enthalpy(const Gas_data &Q, int isp, int itm)
{
    return TBM_->eval_modal_enthalpy_isp(Q, EOS_, isp, itm);
}

double
Composite_gas_model::
s_modal_Cv(Gas_data &Q, int itm)
{
    return TBM_->eval_modal_Cv(Q, EOS_, itm);
}

void
Composite_gas_model::
initialise_ideal_gas(lua_State *L)
{
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Composite_gas_model::initialise_ideal_gas():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    int nsp = lua_objlen(L, -1);
    set_number_of_species(nsp);

    // For an "ideal" gas, nmodes = 1
    set_number_of_modes(1);

    //    cout << "initialising equation of state...\n";
    EOS_ = (Equation_of_state*) new Perfect_gas(L);

    //cout << "initialising thermal behaviour model...\n";
    TBM_ = (Thermal_behaviour_model*) new Constant_specific_heats(L);

    //cout << "initialising transport coefficients model...\n";
    TCM_ = (Transport_coefficients_model*) new Wilke_mixing_rule(L);
    
    if ( nsp > 1 ) {
	cout << "FIX ME: diffusion coefficients model when nsp > 1.\n";
	DCM_ = 0;
	//	DCM_ = (Diffusion_coefficients_model*) new Hard_sphere_diffusion_coefficients(L);
	exit(1);
    }
    else
	DCM_ = new No_diffusion_coefficients();
    
    //cout << "initialising sound speed model...\n";
    SSM_ = create_sound_speed_model( "equilibrium", *this, L );
}

Gas_model* create_composite_gas_model(string cfile)
{
    return new Composite_gas_model(cfile);
}
