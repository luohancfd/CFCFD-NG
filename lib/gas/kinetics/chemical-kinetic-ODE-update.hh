// Author: Rowan J. Gollan
// Date: 12-Sep-2008

#ifndef CHEMICAL_KINETIC_ODE_UPDATE_HH
#define CHEMICAL_KINETIC_ODE_UPDATE_HH

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
#include "reaction-update.hh"
#include "chemical-kinetic-system.hh"
#include "reaction.hh"

class Chemical_kinetic_ODE_update : public Reaction_update {
public:
    Chemical_kinetic_ODE_update(lua_State *L, Gas_model &g);
    ~Chemical_kinetic_ODE_update();
    
    Chemical_kinetic_system * get_cks_pointer()
    { return cks_; }

private:
    int s_update_state(Gas_data &Q, double dt, double &dt_suggest, Gas_model *gm=0);
    int s_rate_of_change(Gas_data &Q, std::vector<double> &dcdt);
    int s_eval_chemistry_energy_coupling_source_terms( Gas_data &Q, std::vector<double> &dedt );

    int perform_increment(Gas_data &Q, double dt, double &dt_suggest);
    int estimate_appropriate_subcycle(double t_interval, double dt_suggest,
				      double &dt_sub, int &no_substeps);
    double T_upper_limit_; // Above this temperature,
                           // the reactions are NOT computed.
    double T_lower_limit_; // Below this temperature,
                           // the reactions are NOT computed.
    OdeSolver *ode_solver_;
    Chemical_kinetic_system *cks_;
    std::valarray<double> yin_, yout_, ydot_;
    std::vector<double> M_;
    Gas_data *Q_save_;
};

Reaction_update* create_Chemical_kinetic_ODE_update(lua_State *L, Gas_model &g);

Chemical_kinetic_system * create_chemical_kinetic_system( std::string cfile, Gas_model &g );

#endif
