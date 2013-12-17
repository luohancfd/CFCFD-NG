// Author: Daniel F. Potter
// Date: 10-May-2010

#ifndef THERMAL_EQUILIBRIUM_MECHANISM_HH
#define THERMAL_EQUILIBRIUM_MECHANISM_HH

#include <string>
// #include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../models/chemical-species.hh"
#include "../models/chemical-species-library.hh"

class Thermal_equilibrium_mechanism {
public:
    Thermal_equilibrium_mechanism( lua_State * L );
     ~Thermal_equilibrium_mechanism();

    void apply( Gas_data &Q )
    { return s_apply( Q ); }

private:
    void s_apply( Gas_data &Q );
    
private:
    int isp_;
    Species_energy_mode * mode_;
    int iT_;
    int iT_eq_;
};

#endif
