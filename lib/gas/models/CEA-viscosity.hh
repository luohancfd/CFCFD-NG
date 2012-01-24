// Author: Rowan J. Gollan
// Date: 05=Nov-2008
// Place: Hampton, Virginia, USA

#ifndef CEA_VISCOSITY_HH
#define CEA_VISCOSITY_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "viscosity-model.hh"
#include "CEA_curves.hh"

class CEA_viscosity : public Viscosity_model {
public:
    CEA_viscosity(lua_State *L);
    ~CEA_viscosity();
private:
    std::vector<CEA_transport_params> curve_;

    double s_eval_viscosity(const Gas_data &Q);
};

#endif
