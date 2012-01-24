// Author: Rowan J. Gollan
// Date: 30-July-2008

#ifndef HARD_SPHERE_DCM_HH
#define HARD_SPHERE_DCM_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "diffusion-coefficients-model.hh"

class Hard_sphere_dcm : public Diffusion_coefficients_model {
public:
    Hard_sphere_dcm(lua_State *L);
    ~Hard_sphere_dcm();
private:
    double R_mix_;           // gas constant for mixture
    std::vector<double> d_;  // hard-sphere diameters of components, m
    std::vector<double> m_;  // masses of single particle, kg
    std::vector<double> R_;  // gas constants for componentes

    int s_eval_diffusion_coefficients(Gas_data &Q);
    double calculate_D_AB(Gas_data &Q, int i, int j);
};

#endif
