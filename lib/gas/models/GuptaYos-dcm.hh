// Author: Daniel F. Potter
// Date: 10-December-2009

#ifndef GUPTAYOS_DCM_HH
#define GUPTAYOS_DCM_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "diffusion-coefficients-model.hh"
#include "binary-interaction.hh"
#include "chemical-species.hh"

class GuptaYos_dcm : public Diffusion_coefficients_model {
public:
    GuptaYos_dcm(lua_State *L);
    ~GuptaYos_dcm();
private:
    std::vector<Binary_interaction*> unique_BIs_;
    std::vector< std::vector<Binary_interaction*> > BI_table_;

    int nsp_;
    std::vector<int> Z_;
    std::vector<double> m_;
    std::vector<double> x_;
    double ignore_mole_fraction_;

    int s_eval_diffusion_coefficients(Gas_data &Q);
};

#endif
