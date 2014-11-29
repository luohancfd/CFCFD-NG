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
    lua_pushnumber(L, cell->iLength); lua_setfield(L, -2, "iLength");
    lua_pushnumber(L, cell->jLength); lua_setfield(L, -2, "jLength");
    lua_pushnumber(L, cell->kLength); lua_setfield(L, -2, "kLength");
    lua_pushnumber(L, cell->volume[0]); lua_setfield(L, -2, "vol");
    lua_pushnumber(L, fs.gas->p); lua_setfield(L, -2, "p");
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

void create_table_for_fs(lua_State *L, FlowState &fs, Gas_model &gmodel)
{
    lua_newtable(L);
    lua_pushnumber(L, fs.gas->p); lua_setfield(L, -2, "p");
    lua_pushnumber(L, fs.gas->rho); lua_setfield(L, -2, "rho"); 
    lua_pushnumber(L, fs.vel.x); lua_setfield(L, -2, "u"); 
    lua_pushnumber(L, fs.vel.y); lua_setfield(L, -2, "v");
    lua_pushnumber(L, fs.vel.z); lua_setfield(L, -2, "w");
    lua_pushnumber(L, fs.gas->a); lua_setfield(L, -2, "a");
    lua_pushnumber(L, fs.gas->mu); lua_setfield(L, -2, "mu");
    lua_pushnumber(L, fs.mu_t); lua_setfield(L, -2, "mu_t");
    lua_pushnumber(L, fs.k_t); lua_setfield(L, -2, "k_t");
    lua_pushnumber(L, fs.tke); lua_setfield(L, -2, "tke");
    lua_pushnumber(L, fs.omega); lua_setfield(L, -2, "omega");
    lua_pushinteger(L, fs.S); lua_setfield(L, -2, "S");
    lua_newtable(L); // A table for the individual temperatures
    size_t nmodes = gmodel.get_number_of_modes();
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	lua_pushinteger(L, static_cast<int>(imode));
	lua_pushnumber(L, fs.gas->T[imode]);
	lua_settable(L, -3);
    }
    // At this point, the table of temperatures should be TOS.
    lua_setfield(L, -2, "T");
    lua_newtable(L); // A table for the individual conductivities
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	lua_pushinteger(L, static_cast<int>(imode));
	lua_pushnumber(L, fs.gas->k[imode]);
	lua_settable(L, -3);
    }
    // At this point, the table of conductivities should be TOS.
    lua_setfield(L, -2, "k");
    lua_newtable(L); // Another table for the mass fractions
    size_t nsp = gmodel.get_number_of_species();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, static_cast<int>(isp));
	lua_pushnumber(L, fs.gas->massf[isp]);
	lua_settable(L, -3);
    }
    // At this point, the table of mass fractions should be TOS.
    lua_setfield(L, -2, "massf");
    return;
}

void create_table_for_iface(lua_State *L, FV_Interface &iface, Gas_model &gmodel)
{
    FlowState &fs = *(iface.fs);
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, iface.pos.x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, iface.pos.y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, iface.pos.z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, iface.length); lua_setfield(L, -2, "length");
    // Flow state information
    lua_pushnumber(L, fs.gas->p); lua_setfield(L, -2, "p");
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
    size_t nmodes = gmodel.get_number_of_modes();
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	lua_pushinteger(L, static_cast<int>(imode));
	lua_pushnumber(L, fs.gas->T[imode]);
	lua_settable(L, -3);
    }
    // At this point, the table of temperatures should be TOS.
    lua_setfield(L, -2, "T");
    lua_newtable(L); // Another table for the mass fractions
    size_t nsp = gmodel.get_number_of_species();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, static_cast<int>(isp));
	lua_pushnumber(L, fs.gas->massf[isp]);
	lua_settable(L, -3);
    }
    // At this point, the table of mass fractions should be TOS.
    lua_setfield(L, -2, "massf");
    return;
}

int luafn_sample_i_face(lua_State *L)
{
    // Get arguments from stack.
    size_t jb = static_cast<size_t>(lua_tointeger(L, 1));
    size_t i = static_cast<size_t>(lua_tointeger(L, 2));
    size_t j = static_cast<size_t>(lua_tointeger(L, 3));
    size_t k = static_cast<size_t>(lua_tointeger(L, 4));

    Gas_model *gmodel = get_gas_model_ptr();
    Block *bdp= get_block_data_ptr(jb);
    FV_Interface *iface = bdp->get_ifi(i,j,k);
    create_table_for_iface(L, *iface, *gmodel);
    return 1; // Just the one table is left on Lua stack
}

int luafn_sample_j_face(lua_State *L)
{
    // Get arguments from stack.
    size_t jb = static_cast<size_t>(lua_tointeger(L, 1));
    size_t i = static_cast<size_t>(lua_tointeger(L, 2));
    size_t j = static_cast<size_t>(lua_tointeger(L, 3));
    size_t k = static_cast<size_t>(lua_tointeger(L, 4));

    Gas_model *gmodel = get_gas_model_ptr();
    Block *bdp= get_block_data_ptr(jb);
    FV_Interface *iface = bdp->get_ifj(i,j,k);
    create_table_for_iface(L, *iface, *gmodel);
    return 1; // Just the one table is left on Lua stack
}

int luafn_sample_k_face(lua_State *L)
{
    // Get arguments from stack.
    size_t jb = static_cast<size_t>(lua_tointeger(L, 1));
    size_t i = static_cast<size_t>(lua_tointeger(L, 2));
    size_t j = static_cast<size_t>(lua_tointeger(L, 3));
    size_t k = static_cast<size_t>(lua_tointeger(L, 4));

    Gas_model *gmodel = get_gas_model_ptr();
    Block *bdp= get_block_data_ptr(jb);
    FV_Interface *iface = bdp->get_ifk(i,j,k);
    create_table_for_iface(L, *iface, *gmodel);
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

    int flag = gmodel->eval_transport_coefficients(Q);
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
    // Assume gas_data is at top of stack and store this index
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

int luafn_molecular_weight(lua_State *L)
{
    Gas_model *gmodel = get_gas_model_ptr();
    // Assume that integer for species is at top of stack
    int isp = luaL_checkint(L, 1);
    double mw = gmodel->molecular_weight(isp);
    // Put mw on top of stack
    lua_pushnumber(L, mw);
    return 1;
}

int luafn_enthalpy(lua_State *L)
{
    // Assume gas_data is a top of stack and store this index
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data Q(gmodel);
    // Gas_data is argument 1. We need to bring it to
    // TOS before using get_table_as_data()
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, Q);
    lua_pop(L, 1);
    // Then expect an integer for species index next
    int isp = luaL_checkinteger(L, 2);
    double h = gmodel->enthalpy(Q, isp);
    // Put h at top of stack
    lua_pushnumber(L, h);
    return 1;
}


int luafn_massf2molef(lua_State *L)
{
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    // Expect a table with mass fractions at top of stack.
    vector<double> massf(nsp, 0.0);
    if ( !lua_istable(L, 1) ) { 
	ostringstream ost;
	ost << "Error in call 'massf2molef()'\n";
	ost << "A table of mass fractions is expected as the first argument.\n";
	luaL_error(L, ost.str().c_str());
    }
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, 1, isp);
	massf[isp] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }
    vector<double> molef(nsp, 0.0);
    convert_massf2molef(massf, gmodel->M(), molef);
    // Push molef into a table.
    lua_newtable(L);
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_pushnumber(L, molef[isp]);
	lua_settable(L, -3);
    }
    return 1;
}

int luafn_massf2conc(lua_State *L)
{
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    // Expect density first
    double rho = luaL_checknumber(L, 1);
    // Now a table of mass fractions
    vector<double> massf(nsp, 0.0);
    if ( !lua_istable(L, 2) ) { 
	ostringstream ost;
	ost << "Error in call 'massf2conc()'\n";
	ost << "A table of mass fractions is expected as the second argument.\n";
	luaL_error(L, ost.str().c_str());
    }
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, 2, isp);
	massf[isp] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }
    vector<double> conc(nsp, 0.0);
    convert_massf2conc(rho, massf, gmodel->M(), conc);
    // Push concentrations into a table.
    lua_newtable(L);
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_pushnumber(L, conc[isp]);
	lua_settable(L, -3);
    }
    return 1;
}

int luafn_molef2massf(lua_State *L)
{
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    // Expect a table with mole fractions at top of stack.
    vector<double> molef(nsp, 0.0);
    if ( !lua_istable(L, 1) ) { 
	ostringstream ost;
	ost << "Error in call 'molef2massf()'\n";
	ost << "A table of mole fractions is expected as the first argument.\n";
	luaL_error(L, ost.str().c_str());
    }
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, 1, isp);
	molef[isp] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }
    vector<double> massf(nsp, 0.0);
    convert_molef2massf(molef, gmodel->M(), massf);
    // Push massf into a table.
    lua_newtable(L);
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_pushnumber(L, massf[isp]);
	lua_settable(L, -3);
    }
    return 1;
}

int luafn_conc2massf(lua_State *L)
{
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    // Expect density first
    double rho = luaL_checknumber(L, 1);
    // Now a table of concentrations fractions
    vector<double> conc(nsp, 0.0);
    if ( !lua_istable(L, 2) ) { 
	ostringstream ost;
	ost << "Error in call 'conc2massf()'\n";
	ost << "A table of concentrations is expected as the second argument.\n";
	luaL_error(L, ost.str().c_str());
    }
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, 2, isp);
	conc[isp] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }
    vector<double> massf(nsp, 0.0);
    convert_conc2massf(rho, conc, gmodel->M(), massf);
    // Push concentrations into a table.
    lua_newtable(L);
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_pushnumber(L, massf[isp]);
	lua_settable(L, -3);
    }
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

int luafn_species_rate_of_change(lua_State *L)
{
    Gas_model *gmodel = get_gas_model_ptr();
    Reaction_update *rupdate = get_reaction_update_ptr();
    if ( rupdate == nullptr ) {
	cout << "ERROR in luafn_species_rate_of_change(): " << endl;
	cout << "There is no reaction scheme intialised in this simulation,\n";
	cout << "so we cannot compute the rates of change for species concentrations.\n";
	cout << "Bailing out!\n";
	exit(UDF_ERROR);
    }
    // Assume gas_data is a top of stack and store this index
    Gas_data Q(gmodel);
    // Expect a gas_data as lua table at top of stack.
    get_table_as_gas_data(L, Q);
    // Initialise a vector to hold concentration rates of change
    vector<double> dcdt(gmodel->get_number_of_species(), 0.0);
    int flag = rupdate->rate_of_change(Q, dcdt);
    if ( flag != SUCCESS ) {
	cout << "luafn_species_rate_of_change(): " << endl;
	cout << "There was a problem calling rate_of_change()." << endl;
	cout << "This was called inside a user-defined function." << endl;
	cout << "Supplied gas_data was:" << endl;
	Q.print_values();
	cout << "Bailing out!" << endl;
	exit(UDF_ERROR);
    }
    // Put table of concentration values are TOS
    push_vector_as_table(L, dcdt);
    return 1;
}

int register_luafns(lua_State *L)
{
    lua_pushcfunction(L, luafn_sample_flow);
    lua_setglobal(L, "sample_flow");
    lua_pushcfunction(L, luafn_sample_i_face);
    lua_setglobal(L, "sample_i_face");
    lua_pushcfunction(L, luafn_sample_j_face);
    lua_setglobal(L, "sample_j_face");
    lua_pushcfunction(L, luafn_sample_k_face);
    lua_setglobal(L, "sample_k_face");
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
    lua_pushcfunction(L, luafn_molecular_weight);
    lua_setglobal(L, "molecular_weight");
    lua_pushcfunction(L, luafn_enthalpy);
    lua_setglobal(L, "enthalpy");
    lua_pushcfunction(L, luafn_massf2molef);
    lua_setglobal(L, "massf2molef");
    lua_pushcfunction(L, luafn_massf2conc);
    lua_setglobal(L, "massf2conc");
    lua_pushcfunction(L, luafn_molef2massf);
    lua_setglobal(L, "molef2massf");
    lua_pushcfunction(L, luafn_conc2massf);
    lua_setglobal(L, "conc2massf");
    lua_pushcfunction(L, luafn_species_rate_of_change);
    lua_setglobal(L, "species_rate_of_change");
    // Set some of the physical constants
    lua_pushnumber(L, PC_R_u);
    lua_setglobal(L, "PC_R_u");
    lua_pushnumber(L, PC_P_atm);
    lua_setglobal(L, "PC_P_atm");
    lua_pushnumber(L, PC_Avogadro);
    lua_setglobal(L, "PC_Avogadro");
    lua_pushnumber(L, PC_k_SI);
    lua_setglobal(L, "PC_k_SI");

    return 0;
}

