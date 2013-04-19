// lua_service_for_e3.cxx
// Author: Rowan J. Gollan
// Date: 13-Mar-2013
// Place: The University of Queensland
//
// History: 13-Mar-2013
//          Refactored some service functions from bc_user_defined.hh
//

#include <iostream>

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "../../../lib/gas/models/lservice_gas_data.hh"
#include "kernel.hh"
#include "block.hh"
#include "cell.hh"
#include "lua_service_for_e3.hh"

using namespace std;

int luafn_sample_flow(lua_State *L)
{
    // Get arguments from stack.
    size_t jb = static_cast<size_t>(lua_tointeger(L, 1));
    size_t i = static_cast<size_t>(lua_tointeger(L, 2));
    size_t j = static_cast<size_t>(lua_tointeger(L, 3));
    size_t k = static_cast<size_t>(lua_tointeger(L, 4));

    Gas_model *gmodel = get_gas_model_ptr();
    Block *bdp= get_block_data_ptr(jb);
    FV_Cell *cell = bdp->get_cell(i,j,k);
    FlowState &fs = *(cell->fs);

    // Return the interesting data in a table with named fields.
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, cell->pos[0].x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, cell->pos[0].y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, cell->pos[0].z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, cell->volume[0]); lua_setfield(L, -2, "vol");
    lua_pushnumber(L, fs.gas->p); lua_setfield(L, -2, "p");
    // lua_pushnumber(L, fs.gas->T[0]); lua_setfield(L, -2, "T_wall"); // FIX-ME: check not needed.
    lua_pushnumber(L, fs.gas->rho); lua_setfield(L, -2, "rho"); 
    lua_pushnumber(L, fs.vel.x); lua_setfield(L, -2, "u"); 
    lua_pushnumber(L, fs.vel.y); lua_setfield(L, -2, "v");
    lua_pushnumber(L, fs.vel.z); lua_setfield(L, -2, "w");
    lua_pushnumber(L, fs.gas->a); lua_setfield(L, -2, "a");
    lua_pushnumber(L, fs.mu_t); lua_setfield(L, -2, "mu_t");
    lua_pushnumber(L, fs.k_t); lua_setfield(L, -2, "k_t");
    lua_pushnumber(L, fs.tke); lua_setfield(L, -2, "tke");
    lua_pushnumber(L, fs.omega); lua_setfield(L, -2, "omega");
    lua_pushinteger(L, fs.S); lua_setfield(L, -2, "S");
    lua_newtable(L); // A table for the individual temperatures
    size_t nmodes = gmodel->get_number_of_modes();
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	lua_pushinteger(L, static_cast<int>(imode));
	lua_pushnumber(L, fs.gas->T[imode]);
	lua_settable(L, -3);
    }
    // At this point, the table of temperatures should be TOS.
    lua_setfield(L, -2, "T");
    lua_newtable(L); // Another table for the mass fractions
    size_t nsp = gmodel->get_number_of_species();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, static_cast<int>(isp));
	lua_pushnumber(L, fs.gas->massf[isp]);
	lua_settable(L, -3);
    }
    // At this point, the table of mass fractions should be TOS.
    lua_setfield(L, -2, "massf");
    // Parameters for computing density and energy residuals
    lua_pushnumber(L, cell->U[0]->mass); lua_setfield(L, -2, "rho_old");
    // Have commented out the following line because I can't guarantee
    // which time-level, ftl, is appropriate for this particular call. 
    // lua_pushnumber(L, cell->U[ftl]->total_energy); lua_setfield(L, -2, "rE");
    lua_pushnumber(L, cell->U[0]->total_energy); lua_setfield(L, -2, "rE_old");
    
    return 1; // Just the one table is left on Lua stack
}

int luafn_locate_cell(lua_State *L)
{
    double x = lua_tonumber(L, 1);
    double y = lua_tonumber(L, 2);
    double z = 0.0;
    global_data *gdp = get_global_data_ptr();
    if ( gdp->dimensions == 3 ) {
	z = lua_tonumber(L, 3);
    }
    size_t jb, i, j, k;
    int found_cell = locate_cell(x, y, z, &jb, &i, &j, &k, 0);
    if ( !found_cell ) {
	found_cell = find_nearest_cell(x, y, z, &jb, &i, &j, &k, 0);
	if ( !found_cell ) {
	    printf("luafn_locate_cell(): no cells near pos=(%g,%g,%g)", x, y, z);
	}
    }
    lua_pushinteger(L, static_cast<int>(jb));
    lua_pushinteger(L, static_cast<int>(i));
    lua_pushinteger(L, static_cast<int>(j));
    lua_pushinteger(L, static_cast<int>(k));
    lua_pushinteger(L, found_cell);
    return 5; // We leave jb, i, j, k and found_cell on the stack.
}

int luafn_create_empty_gas_table(lua_State *L)
{
    Gas_model* gmodel = get_gas_model_ptr();
    return create_empty_gas_table(L, *gmodel);
}

int luafn_eval_thermo_state_pT(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    int tindex = lua_gettop(L);
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    get_table_as_gas_data(L, Q);
    int flag = apply_gas_method(&Gas_model::eval_thermo_state_pT, Q);
    if ( flag != SUCCESS ) {
	cout << "luafn_eval_thermo_state_pT(): " << endl;
	cout << "There was a problem calling eval_thermo_state_pT()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Modify table at top of stack
    set_gas_data_at_table(L, tindex, Q);
    return 0;
}

int luafn_eval_thermo_state_rhoe(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    int tindex = lua_gettop(L);
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int flag = apply_gas_method(&Gas_model::eval_thermo_state_rhoe, Q);
    if ( flag != SUCCESS ) {
	cout << "luafn_eval_thermo_state_rhoe(): " << endl;
	cout << "There was a problem calling eval_thermo_state_rhoe()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Modify table at top of stack
    set_gas_data_at_table(L, tindex, Q);
    return 0;
}

int luafn_eval_thermo_state_rhoT(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    int tindex = lua_gettop(L);
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int flag = apply_gas_method(&Gas_model::eval_thermo_state_rhoT, Q);
    if ( flag != SUCCESS ) {
	cout << "luafn_eval_thermo_state_rhoT(): " << endl;
	cout << "There was a problem calling eval_thermo_state_rhoT()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Modify table at top of stack
    set_gas_data_at_table(L, tindex, Q);
    return 0;
}

int luafn_eval_thermo_state_rhop(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    int tindex = lua_gettop(L);
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int flag = apply_gas_method(&Gas_model::eval_thermo_state_rhop, Q);
    if ( flag != SUCCESS ) {
	cout << "luafn_eval_thermo_state_rhop(): " << endl;
	cout << "There was a problem calling eval_thermo_state_rhop()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Modify table at top of stack
    set_gas_data_at_table(L, tindex, Q);
    return 0;
}

int luafn_eval_sound_speed(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    int tindex = lua_gettop(L);
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int flag = apply_gas_method(&Gas_model::eval_sound_speed, Q);
    if ( flag != SUCCESS ) {
	cout << "luafn_eval_thermo_sound_speed(): " << endl;
	cout << "There was a problem calling eval_sound_speed()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Modify table at top of stack
    set_gas_data_at_table(L, tindex, Q);
    return 0;
}

int luafn_eval_transport_coefficients(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    int tindex = lua_gettop(L);
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int flag = apply_gas_method(&Gas_model::eval_transport_coefficients, Q);
    if ( flag != SUCCESS ) {
	cout << "luafn_eval_transport_coefficients(): " << endl;
	cout << "There was a problem calling eval_transport_coefficients()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Modify table at top of stack
    set_gas_data_at_table(L, tindex, Q);
    return 0;
}

int luafn_eval_diffusion_coefficients(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    int tindex = lua_gettop(L);
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int flag = apply_gas_method(&Gas_model::eval_diffusion_coefficients, Q);
    if ( flag != SUCCESS ) {
	cout << "luafn_eval_diffusion_coefficients(): " << endl;
	cout << "There was a problem calling eval_diffusion_coefficients()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Modify table at top of stack
    set_gas_data_at_table(L, tindex, Q);
    return 0;
}

int luafn_eval_Cv(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int status;
    double Cv = apply_gas_method(&Gas_model::Cv, Q, status);
    if ( status != SUCCESS ) {
	cout << "luafn_eval_Cv(): " << endl;
	cout << "There was a problem calling eval_Cv()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Put Cv at top of stack
    lua_pushnumber(L, Cv);
    return 1;
}

int luafn_eval_Cp(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int status;
    double Cp = apply_gas_method(&Gas_model::Cp, Q, status);
    if ( status != SUCCESS ) {
	cout << "luafn_eval_Cp(): " << endl;
	cout << "There was a problem calling eval_Cp()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Put Cp at top of stack
    lua_pushnumber(L, Cp);
    return 1;
}

int luafn_eval_R(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int status;
    double R = apply_gas_method(&Gas_model::R, Q, status);
    if ( status != SUCCESS ) {
	cout << "luafn_eval_R(): " << endl;
	cout << "There was a problem calling eval_R()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Put R at top of stack
    lua_pushnumber(L, R);
    return 1;
}

int luafn_eval_gamma(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    int status;
    double gamma = apply_gas_method(&Gas_model::gamma, Q, status);
    if ( status != SUCCESS ) {
	cout << "luafn_eval_R(): " << endl;
	cout << "There was a problem calling eval_gamma()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Put gamma at top of stack
    lua_pushnumber(L, gamma);
    return 1;
}

int apply_gas_method(Gas_model_Method_gas_data f, Gas_data &Q)
{
    Gas_model *gmodel = get_gas_model_ptr();
    return CALL_MEMBER_FN(*gmodel,f)(Q);
}

double apply_gas_method(Gas_model_Method_gas_data_int f, Gas_data &Q, int &status)
{
    Gas_model *gmodel = get_gas_model_ptr();
    return CALL_MEMBER_FN(*gmodel,f)(Q, status);
}

int register_luafns(lua_State *L)
{
    lua_pushcfunction(L, luafn_sample_flow);
    lua_setglobal(L, "sample_flow");
    lua_pushcfunction(L, luafn_locate_cell);
    lua_setglobal(L, "locate_cell");
    lua_pushcfunction(L, luafn_create_empty_gas_table);
    lua_setglobal(L, "create_empty_gas_table");
    lua_pushcfunction(L, luafn_eval_thermo_state_pT);
    lua_setglobal(L, "eval_thermo_state_pT");
    lua_pushcfunction(L, luafn_eval_thermo_state_rhoe);
    lua_setglobal(L, "eval_thermo_state_rhoe");
    lua_pushcfunction(L, luafn_eval_thermo_state_rhoT);
    lua_setglobal(L, "eval_thermo_state_rhoT");
    lua_pushcfunction(L, luafn_eval_thermo_state_rhop);
    lua_setglobal(L, "eval_thermo_state_rhop");
    lua_pushcfunction(L, luafn_eval_sound_speed);
    lua_setglobal(L, "eval_sound_speed");
    lua_pushcfunction(L, luafn_eval_transport_coefficients);
    lua_setglobal(L, "eval_transport_coefficients");
    lua_pushcfunction(L, luafn_eval_diffusion_coefficients);
    lua_setglobal(L, "eval_diffusion_coefficients");
    lua_pushcfunction(L, luafn_eval_Cv);
    lua_setglobal(L, "eval_Cv");
    lua_pushcfunction(L, luafn_eval_Cp);
    lua_setglobal(L, "eval_Cp");
    lua_pushcfunction(L, luafn_eval_R);
    lua_setglobal(L, "eval_R");
    lua_pushcfunction(L, luafn_eval_gamma);
    lua_setglobal(L, "eval_gamma");
    // Set some of the physical constants
    lua_pushnumber(L, PC_R_u);
    lua_setglobal(L, "PC_R_u");
    lua_pushnumber(L, PC_P_atm);
    lua_setglobal(L, "PC_P_atm");

    return 0;
}

