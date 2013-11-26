// Author: Daniel F. Potter
// Date: 18-Nov-2009

#ifndef ENERGY_EXCHANGE_ODE_UPDATE_HH
#define ENERGY_EXCHANGE_ODE_UPDATE_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../models/gas-model.hh"
#include "../../nm/source/ode_solver.hh"
#include "../../nm/source/ode_system.hh"

#include "energy-exchange-update.hh"
#include "energy-exchange-system.hh"
#include "thermal-equilibrium-mechanism.hh"

class Energy_exchange_ODE_update : public Energy_exchange_update {
public:
    Energy_exchange_ODE_update(lua_State *L, Gas_model &g);
    ~Energy_exchange_ODE_update();
    
    Energy_exchange_system * get_ees_pointer()
    { return ees_; }

private:
    int s_update_state( Gas_data &Q, double dt, double &dt_suggest, Gas_model *gm);
    int s_rate_of_change(Gas_data &Q, std::vector<double> &dedt);

    int perform_increment(Gas_data &Q, double dt, double &dt_suggest);
    int estimate_appropriate_subcycle(double t_interval, double dt_suggest,
				      double &dt_sub, int &no_substeps);
    
    double T_upper_limit_; // Above this temperature,
                           // energy exchange is NOT computed.
    double T_lower_limit_; // Below this temperature,
                           // energy exchange is NOT computed.
    OdeSolver *ode_solver_;
    Energy_exchange_system *ees_;
    std::vector<double> yin_, yout_, ydot_;
    Gas_data *Q_save_;
    std::vector<double> molef_;
    Gas_model * g_;
    std::vector<Thermal_equilibrium_mechanism*> eq_mechs_;
};

Energy_exchange_update* create_Energy_exchange_ODE_update(lua_State *L, Gas_model &g);

#endif
