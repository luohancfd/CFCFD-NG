// Author: Daniel F. Potter
// Date: 18-Nov-2009

#ifndef ENERGY_EXCHANGE_RATE_HH
#define ENERGY_EXCHANGE_RATE_HH

#include <vector>
#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"

#include "energy-exchange-mechanism.hh"

class Energy_exchange_rate {
public:
    Energy_exchange_rate(lua_State *L);

    ~Energy_exchange_rate();

    int compute_all_relaxation_times(Gas_data &Q, std::vector<double> &molef);

    double compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);

    bool is_equilibriated()
    { return equilibriated_; }

    int equilibriated_with_mode()
    { return imode_eq_; }

private:
    bool equilibriated_;
    int imode_eq_;
    std::vector<Energy_exchange_mechanism*> ee_mech_;
};

Energy_exchange_rate * create_energy_exchange_rate(lua_State *L);

#endif
