// Author: Rowan J. Gollan
// Date: 09-Jul-2008

#ifndef COMPOSITE_GAS_MODEL_HH
#define COMPOSITE_GAS_MODEL_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "gas-model.hh"
#include "equation-of-state.hh"
#include "thermal-behaviour-model.hh"
#include "transport-coefficients-model.hh"
#include "diffusion-coefficients-model.hh"
#include "sound-speed-model.hh"

class Composite_gas_model : public Gas_model {
public:
    Composite_gas_model(std::string cfile);
    ~Composite_gas_model();

private:
    Equation_of_state *EOS_;
    Thermal_behaviour_model *TBM_;
    Transport_coefficients_model *TCM_;
    Diffusion_coefficients_model *DCM_;
    Sound_speed_model *SSM_;

    int s_decode_conserved_energy(Gas_data &Q, const std::vector<double> &rhoe);
    int s_encode_conserved_energy(const Gas_data &Q, std::vector<double> &rhoe);
    int s_eval_thermo_state_rhoe(Gas_data &Q);
    int s_eval_thermo_state_pT(Gas_data &Q);
    int s_eval_thermo_state_rhoT(Gas_data &Q);
    int s_eval_thermo_state_rhop(Gas_data &Q);
    int s_eval_sound_speed(Gas_data &Q);
    int s_eval_transport_coefficients(Gas_data &Q);
    int s_eval_diffusion_coefficients(Gas_data &Q);
    double s_dTdp_const_rho(const Gas_data &Q, int &status);
    double s_dTdrho_const_p(const Gas_data &Q, int &status);
    double s_dpdrho_const_T(const Gas_data &Q, int &status);
    double s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status);
    double s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status);
    double s_dedT_const_v(const Gas_data &Q, int &status);
    double s_dhdT_const_p(const Gas_data &Q, int &status);
    double s_internal_energy(const Gas_data &Q, int isp);
    double s_enthalpy(const Gas_data &Q, int isp);
    double s_entropy(const Gas_data &Q, int isp);
    double s_modal_enthalpy( double T, int isp, int itm);
    double s_modal_Cv(Gas_data &Q, int itm);

    void initialise_ideal_gas(lua_State *L);
};

Gas_model* create_composite_gas_model(std::string cfile);

#endif
