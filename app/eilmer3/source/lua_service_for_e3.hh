// lua_service_for_e3.hh
// Author: Rowan J. Gollan
// Date: 13-Mar-2013
// Place: The University of Queensland
//
// History: 13-Mar-2013
//          Refactored some service functions from bc_user_defined.hh
//

#ifndef LUA_SERVICE_FOR_E3_HEADER
#define LUA_SERVICE_FOR_E3_HEADER

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

typedef int (Gas_model::*Gas_model_Method_gas_data)(Gas_data &);
typedef double (Gas_model::*Gas_model_Method_gas_data_int)(const Gas_data &, int &);
#define CALL_MEMBER_FN(object,ptrToMember) ((object).*(ptrToMember))

int luafn_sample_flow(lua_State *L);
void create_table_for_fs(lua_State *L, FlowState &fs, Gas_model &gmodel);
void create_table_for_iface(lua_State *L, FV_Interface &iface, Gas_model &gmodel);
int luafn_sample_i_face(lua_State *L);
int luafn_sample_j_face(lua_State *L);
int luafn_sample_k_face(lua_State *L);
int luafn_locate_cell(lua_State *L);
int luafn_create_empty_gas_table(lua_State *L);
int luafn_eval_thermo_state_pT(lua_State *L);
int luafn_eval_thermo_state_rhoe(lua_State *L);
int luafn_eval_thermo_state_rhoT(lua_State *L);
int luafn_eval_thermo_state_rhop(lua_State *L);
int luafn_eval_sound_speed(lua_State *L);
int luafn_eval_transport_coefficients(lua_State *L);
int luafn_eval_diffusion_coefficients(lua_State *L);
int luafn_eval_Cv(lua_State *L);
int luafn_eval_Cp(lua_State *L);
int luafn_eval_R(lua_State *L);
int luafn_eval_gamma(lua_State *L);
int luafn_molecular_weight(lua_State *L);
int luafn_massf2molef(lua_State *L);
int luafn_massf2conc(lua_State *L);
int luafn_molef2massf(lua_State *L);
int luafn_conc2massf(lua_State *L);
int apply_gas_method(Gas_model_Method_gas_data f, Gas_data &Q);
double apply_gas_method(Gas_model_Method_gas_data_int f, Gas_data &Q, int &status);
int luafn_species_rate_of_change(lua_State *L);
int register_luafns(lua_State *L);

#endif
