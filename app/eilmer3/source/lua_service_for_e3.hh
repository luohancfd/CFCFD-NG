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
#define CALL_MEMBER_FN(object,ptrToMember) ((object).*(ptrToMember))

int luafn_sample_flow(lua_State *L);
int luafn_locate_cell(lua_State *L);
int luafn_create_empty_gas_table(lua_State *L);
int luafn_eval_thermo_state_pT(lua_State *L);
int luafn_eval_thermo_state_rhoe(lua_State *L);
int luafn_eval_thermo_state_rhoT(lua_State *L);
int luafn_eval_thermo_state_rhop(lua_State *L);
int apply_gas_method(Gas_model_Method_gas_data f, Gas_data &Q);
int register_luafns(lua_State *L);

#endif
