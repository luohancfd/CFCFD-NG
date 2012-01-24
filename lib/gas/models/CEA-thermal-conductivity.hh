// Author: Rowan J. Gollan
// Date: 30-July-2008

#ifndef CEA_THERMAL_CONDUCTIVITY_HH
#define CEA_THERMAL_CONDUCTIVITY_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "thermal-conductivity-model.hh"
#include "CEA_curves.hh"

class CEA_thermal_conductivity : public Thermal_conductivity_model {
public:
    CEA_thermal_conductivity(lua_State *L);
    ~CEA_thermal_conductivity();
private:
    std::vector<CEA_transport_params> curve_;

    double s_eval_thermal_conductivity(const Gas_data &Q);
};

#endif
