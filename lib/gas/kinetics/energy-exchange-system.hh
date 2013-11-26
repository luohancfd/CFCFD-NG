// Author: Daniel F. Potter
// Date: 18-Nov-2009

#ifndef ENERGY_EXCHANGE_SYSTEM_HH
#define ENERGY_EXCHANGE_SYSTEM_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/ode_system.hh"
#include "../models/gas-model.hh"
#include "../models/gas_data.hh"

#include "energy-exchange-rate.hh"

class Energy_exchange_system : public OdeSystem {
public:
    Energy_exchange_system(lua_State *L, Gas_model &g, double error_tol);
    ~Energy_exchange_system();
    
    int eval( const std::vector<double> &y, std::vector<double> &ydot );
    double stepsize_select(const std::vector<double> &y);

    void set_gas_data_ptr(Gas_data &Q)
    { Q_ = &Q; }
    
    void set_molef_ptr(std::vector<double> &molef)
    { molef_ = &molef; }

private:
    Gas_model * g_;
    std::vector<Energy_exchange_rate*> ee_rate_;
    double err_tol_;
    Gas_data *Q_;
    std::vector<double> *molef_;
    std::vector<double> ydot_;
};

#endif
