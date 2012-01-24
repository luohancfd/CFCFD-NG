// Author: Rowan J. Gollan
// Date: 24-Mar-2009
// Place: Poquoson, Virginia, USA
//

#ifndef LGAS_MODEL_HH
#define LGAS_MODEL_HH

#include <string>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/lunar.hh"
#include "gas-model.hh"

class luaGas_model {
public:
    static const char className[];
    static Lunar<luaGas_model>::RegType member_data[];
    static Lunar<luaGas_model>::RegType methods[];
    static Lunar<luaGas_model>::MetaType metamethods[];

    std::string str() const;

    Gas_model* r_pointer()
    { return g_; }

    luaGas_model(lua_State *L);
    ~luaGas_model();

    int get_number_of_species(lua_State *L);
    int get_number_of_modes(lua_State *L);

    int eval_thermo_state_pT(lua_State *L);
    int eval_thermo_state_rhoe(lua_State *L);
    int eval_thermo_state_rhoT(lua_State *L);
    int eval_thermo_state_rhop(lua_State *L);
    
    int eval_sound_speed(lua_State *L);
    int eval_transport_coefficients(lua_State *L);
    int eval_diffusion_coefficients(lua_State *L);
    
    int dTdp_const_rho(lua_State *L);
    int dTdrho_const_p(lua_State *L);
    int dpdrho_const_T(lua_State *L);
    int dedT_const_v(lua_State *L);
    int dhdT_const_p(lua_State *L);

    int Cv(lua_State *L);
    int Cp(lua_State *L);
    int R(lua_State *L);
    int gamma(lua_State *L);
    int Prandtl(lua_State *L);
    int mixture_molecular_weight(lua_State *L);
    int molecular_weight(lua_State *L);
    int internal_energy(lua_State *L);
    int enthalpy(lua_State *L);
    int entropy(lua_State *L);
    int Gibbs_free_energy(lua_State *L);
    int species_name(lua_State *L);

private:
    Gas_model *g_;
    Gas_data *Q_;
};

extern "C" int luaopen_gas(lua_State *L);

#endif
