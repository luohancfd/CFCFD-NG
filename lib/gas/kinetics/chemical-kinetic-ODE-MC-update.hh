// Author: Rowan J. Gollan
// Date: 12-Sep-2008

#ifndef CHEMICAL_KINETIC_ODE_MC_UPDATE_HH
#define CHEMICAL_KINETIC_ODE_MC_UPDATE_HH

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
#include "chemical-kinetic-MC-system.hh"
#include "reaction.hh"

class Chemical_kinetic_ODE_MC_update : public Reaction_update {
public:
    Chemical_kinetic_ODE_MC_update(lua_State *L, Gas_model &g);
    ~Chemical_kinetic_ODE_MC_update();
    
    Chemical_kinetic_MC_system * get_cks_pointer()
    { return cks_; }

private:
    int s_update_state(Gas_data &Q, double dt, double &dt_suggest, Gas_model *gm=0);
    int s_rate_of_change(Gas_data &Q, std::vector<double> &dcdt);
    int s_eval_chemistry_energy_coupling_source_terms( Gas_data &Q, std::vector<double> &dedt );
    int s_get_directional_rates( std::vector<double> &w_f, std::vector<double> &w_b )
    { return cks_->get_directional_rates(w_f,w_b); }

    int perform_increment(Gas_data &Q, double dt, double &dt_suggest);
    int estimate_appropriate_subcycle(double t_interval, double dt_suggest,
				      double &dt_sub, int &no_substeps);
    OdeSolver *ode_solver_;
    Chemical_kinetic_MC_system *cks_;
    std::valarray<double> yin_, yout_, ydot_;
    std::valarray<double> cdot_;
    std::vector<double> c_;
    std::vector<double> M_;
    Gas_data *Q_save_;
    Gas_model *gm_;
};

Reaction_update* create_Chemical_kinetic_ODE_MC_update(lua_State *L, Gas_model &g);

Chemical_kinetic_MC_system * create_chemical_kinetic_MC_system( std::string cfile, Gas_model &g );

#endif
