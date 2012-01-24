// Author: Rowan J. Gollan
// Date: 24-Mar-2009
// Place: Poquoson, Virginia, USA

#include <sstream>
#include <string>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "gas-model.hh"
#include "lgas-model.hh"
#include "lservice_gas_data.hh"
#include "physical_constants.hh"
#include "../kinetics/lreaction-update.hh"
#include "../kinetics/lreaction-rate-coeff.hh"

using namespace std;

luaGas_model::
luaGas_model(lua_State *L)
{
    string cfile(luaL_checkstring(L, 1));
    g_ = create_gas_model(cfile);
    Q_ = new Gas_data(g_);

    // Set up a global species table..
    lua_newtable(L);
    for ( int isp = 0; isp < g_->get_number_of_species(); ++isp ) {
	lua_pushstring(L, g_->species_name(isp).c_str());
	lua_rawseti(L, -2, isp+1);
	lua_pushinteger(L, isp+1);
	lua_setfield(L, -2, g_->species_name(isp).c_str());
    }
    lua_setglobal(L, "species");

}

luaGas_model::
~luaGas_model()
{
    delete g_;
}

const char luaGas_model::className[] = "Gas_model";

#define member_data(class, name) {#name, &class::name}

Lunar<luaGas_model>::RegType luaGas_model::member_data[] = {
    {0, 0}
};

#define method(class, name) {#name, &class::name}

Lunar<luaGas_model>::RegType luaGas_model::methods[] = {
    method(luaGas_model, get_number_of_species),
    method(luaGas_model, get_number_of_modes),
    method(luaGas_model, eval_thermo_state_pT),
    method(luaGas_model, eval_thermo_state_rhoe),
    method(luaGas_model, eval_thermo_state_rhoT),
    method(luaGas_model, eval_thermo_state_rhop),
    method(luaGas_model, eval_sound_speed),
    method(luaGas_model, eval_transport_coefficients),
    method(luaGas_model, eval_diffusion_coefficients),
    method(luaGas_model, dTdp_const_rho),
    method(luaGas_model, dTdrho_const_p),
    method(luaGas_model, dpdrho_const_T),
    method(luaGas_model, dedT_const_v),
    method(luaGas_model, dhdT_const_p),
    method(luaGas_model, Cv),
    method(luaGas_model, Cp),
    method(luaGas_model, R),
    method(luaGas_model, gamma),
    method(luaGas_model, Prandtl),
    method(luaGas_model, mixture_molecular_weight),
    method(luaGas_model, molecular_weight),
    //    method(luaGas_model, internal_energy),
    method(luaGas_model, enthalpy),
    method(luaGas_model, entropy),
    method(luaGas_model, Gibbs_free_energy),
    method(luaGas_model, species_name),
    {0, 0}
};

Lunar<luaGas_model>::MetaType luaGas_model::metamethods[] = {
    {0, 0}
};

string
luaGas_model::
str() const
{
    return string("");
}


int
luaGas_model::
get_number_of_species(lua_State *L)
{
    lua_pushnumber(L, g_->get_number_of_species());
    return 1;
}

int
luaGas_model::
get_number_of_modes(lua_State *L)
{
    lua_pushnumber(L, g_->get_number_of_modes());
    return 1;
}

int
luaGas_model::
eval_thermo_state_pT(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to eval_thermo_state_pT()\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    if ( g_->eval_thermo_state_pT(*Q_) == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    set_gas_data_at_table(L, 1, *Q_);

    return 1;
}

int
luaGas_model::
eval_thermo_state_rhoe(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to eval_thermo_state_rhoe():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    if ( g_->eval_thermo_state_rhoe(*Q_) == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    set_gas_data_at_table(L, 1, *Q_);

    return 1;
}

int
luaGas_model::
eval_thermo_state_rhoT(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to eval_thermo_state_rhoT()\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    if ( g_->eval_thermo_state_rhoT(*Q_) == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    set_gas_data_at_table(L, 1, *Q_);

    return 1;
}

int
luaGas_model::
eval_thermo_state_rhop(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to eval_thermo_state_rhop()\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    if ( g_->eval_thermo_state_rhop(*Q_) == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    set_gas_data_at_table(L, 1, *Q_);

    return 1;
}

int
luaGas_model::
eval_sound_speed(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to eval_sound_speed():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    if ( g_->eval_sound_speed(*Q_) == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    set_gas_data_at_table(L, 1, *Q_);

    return 1;
}

int
luaGas_model::
eval_transport_coefficients(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to eval_transport_coefficients():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    if ( g_->eval_transport_coefficients(*Q_) == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    set_gas_data_at_table(L, 1, *Q_);

    return 1;
}

int
luaGas_model::
eval_diffusion_coefficients(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to eval_diffusion_coefficients():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    if ( g_->eval_diffusion_coefficients(*Q_) == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    set_gas_data_at_table(L, 1, *Q_);

    return 1;
}

int
luaGas_model::
dTdp_const_rho(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to dTdp_const_rho():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double deriv;
    deriv = g_->dTdp_const_rho(*Q_, status);
    lua_pushnumber(L, deriv);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
dTdrho_const_p(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to dTdrho_const_p():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double deriv;
    deriv = g_->dTdrho_const_p(*Q_, status);
    lua_pushnumber(L, deriv);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
dpdrho_const_T(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to dpdrho_const_T():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double deriv;
    deriv = g_->dpdrho_const_T(*Q_, status);
    lua_pushnumber(L, deriv);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
dedT_const_v(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to dedT_const_v():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double deriv;
    deriv = g_->dedT_const_v(*Q_, status);
    lua_pushnumber(L, deriv);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
dhdT_const_p(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to dhdT_const_p():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double deriv;
    deriv = g_->dhdT_const_p(*Q_, status);
    lua_pushnumber(L, deriv);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
Cv(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to Cv():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double val;
    val = g_->Cv(*Q_, status);
    lua_pushnumber(L, val);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
Cp(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to Cp():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double val;
    val = g_->Cp(*Q_, status);
    lua_pushnumber(L, val);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
R(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to R():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double val;
    val = g_->R(*Q_, status);
    lua_pushnumber(L, val);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
gamma(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to gamma():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int status;
    double val;
    val = g_->gamma(*Q_, status);
    lua_pushnumber(L, val);
    if ( status == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaGas_model::
Prandtl(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 3 ) {
	ostringstream ost;
	ost << "Error in call to Prandtl():\n";
	ost << "3 arguments expected: mu, Cp and k.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    double mu = luaL_checknumber(L, 1);
    double Cp = luaL_checknumber(L, 2);
    double k = luaL_checknumber(L, 3);
    double Prandtl = g_->Prandtl(mu, Cp, k);

    lua_pushnumber(L, Prandtl);
    return 1;
}

int
luaGas_model::
mixture_molecular_weight(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to mixture_molecular_weight():\n";
	ost << "1 argument expected, a gas_data table.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);
    
    double MW = g_->mixture_molecular_weight(*Q_);
    lua_pushnumber(L, MW);
    return 1;
}

int
luaGas_model::
molecular_weight(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to molecular_weight():\n";
	ost << "1 argument expected, a species index.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    int isp = luaL_checkint(L, 1);
    
    double MW = g_->molecular_weight(isp-1); // -1 for C++ 0-offset
    lua_pushnumber(L, MW);
    return 1;
}

// int
// luaGas_model::
// internal_energy(lua_State *L)
// {
//     int narg = lua_gettop(L);
//     if ( narg != 2 ) {
// 	ostringstream ost;
// 	ost << "Error in call to internal_energy():\n";
// 	ost << "2 arguments expected, a gas_data table and a species index.\n";
// 	ost << narg << " argument(s) received.\n";
// 	input_error(ost);
//     }

//     lua_pushvalue(L, 1);
//     get_table_as_gas_data(L, Q_);
//     lua_pop(L, 1);
    
//     int isp = luaL_checkint(L, 2);
    
//     double e = g_->internal_energy(Q_, isp-1); // -1 for C++ 0-offset
//     lua_pushnumber(L, e);
//     return 1;
// }

int
luaGas_model::
enthalpy(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 2 ) {
	ostringstream ost;
	ost << "Error in call to enthalpy():\n";
	ost << "2 arguments expected, a gas_data table and a species index.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);
    
    int isp = luaL_checkint(L, 2);
    
    double h = g_->enthalpy(*Q_, isp-1); // -1 for C++ 0-offset
    lua_pushnumber(L, h);
    return 1;
}

int
luaGas_model::
entropy(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 2 ) {
	ostringstream ost;
	ost << "Error in call to entropy():\n";
	ost << "2 arguments expected, a gas_data table and a species index.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);
    
    int isp = luaL_checkint(L, 2);
    
    double s = g_->entropy(*Q_, isp-1); // -1 for C++ 0-offset
    lua_pushnumber(L, s);
    return 1;
}

int
luaGas_model::
Gibbs_free_energy(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 2 ) {
	ostringstream ost;
	ost << "Error in call to Gibbs_free_energy():\n";
	ost << "2 arguments expected, a gas_data table and a species index.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);
    
    int isp = luaL_checkint(L, 2);
    
    double g = g_->Gibbs_free_energy(*Q_, isp-1); // -1 for C++ 0-offset
    lua_pushnumber(L, g);
    return 1;
}

int
luaGas_model::
species_name(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to species_name():\n";
	ost << "1 argument expected, a species index.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }
    
    int isp = luaL_checkint(L, 1);

    string name(g_->species_name(isp-1)); // -1 for C++ 0-offset
    
    lua_pushstring(L, name.c_str());
    return 1;
}

int luaopen_gas(lua_State *L)
{
    lua_newtable(L);
    lua_setglobal(L, "gas");
    
    lua_getglobal(L, "gas");
    int table = lua_gettop(L);
    Lunar<luaGas_model>::Register(L, table);

    lua_pushcfunction(L, create_empty_gas_table);
    lua_setfield(L, table, "gas_data");

    // Set some (useful) physical constants
    lua_pushnumber(L, PC_R_u);
    lua_setfield(L, table, "R_u");
    lua_pushnumber(L, PC_Avogadro);
    lua_setfield(L, table, "Avogadro");
    lua_pushnumber(L, PC_k_SI);
    lua_setfield(L, table, "kB");
    lua_pushnumber(L, PC_T_ref);
    lua_setfield(L, table, "T_ref");
    lua_pushnumber(L, PC_P_atm);
    lua_setfield(L, table, "P_atm");

    open_reaction_update(L, table);
    open_reaction_rate_coefficient(L, table);
	
    return 1;
}
