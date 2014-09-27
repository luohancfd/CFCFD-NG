// Author: Rowan J. Gollan
// Date: 03-Jul-2008

#ifndef UD_GAS_MODEL_HH
#define UD_GAS_MODEL_HH

#include <string>
#include <vector>
#include <sstream>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "gas-model.hh"

/** \brief A gas model where the behaviour is supplied by a user-supplied Lua file.
 *
 * This class gives the user some flexibility in implementing some specialised
 * gas model behaviour for some one-off situations.  If the user's model is
 * widely useful, it might be a candidate to be coded as a C++ class. In this
 * sense, the Lua script implementation can be used for rapid prototyping of
 * a new gas model implementation.
 **/

class UD_gas_model : public Gas_model {
public:
    UD_gas_model(std::string cfile);
    ~UD_gas_model();

private:
    lua_State *L_;
    
    int s_eval_thermo_state_rhoe(Gas_data &Q);
    int s_eval_thermo_state_pT(Gas_data &Q);
    int s_eval_thermo_state_rhoT(Gas_data &Q);
    int s_eval_thermo_state_rhop(Gas_data &Q);
    int s_eval_transport_coefficients(Gas_data &Q, Gas_model *gmodel);
    int s_eval_diffusion_coefficients(Gas_data &Q);
    double s_dTdp_const_rho(const Gas_data &Q, int &status);
    double s_dTdrho_const_p(const Gas_data &Q, int &status);
    double s_dpdrho_const_T(const Gas_data &Q, int &status);
    double s_dedT_const_v(const Gas_data &Q, int &status);
    double s_dhdT_const_p(const Gas_data &Q, int &status);
    double s_gas_constant(const Gas_data &Q, int &status);
    double s_molecular_weight(int isp);
    double s_internal_energy(const Gas_data &Q, int isp);
    double s_enthalpy(const Gas_data &Q, int isp);
    double s_entropy(const Gas_data &Q, int isp);
};

Gas_model* create_user_defined_gas_model(std::string cfile);


#endif
