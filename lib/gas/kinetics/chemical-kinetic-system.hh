// Author: Rowan J. Gollan
// Date: 12-Sep-2008

#ifndef CHEMICAL_KINETIC_SYSTEM_HH
#define CHEMICAL_KINETIC_SYSTEM_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/ode_system.hh"
#include "../models/gas_data.hh"
#include "../models/gas-model.hh"
#include "reaction.hh"

class Chemical_kinetic_system : public OdeSystem {
public:
    Chemical_kinetic_system(lua_State *L, Gas_model &g, double error_tol, double T_upper, double T_lower);
    Chemical_kinetic_system(std::string cfile, Gas_model &g);

    ~Chemical_kinetic_system();
    
    int eval(const std::vector<double> &y, std::vector<double> &ydot);
    int eval_split(const std::vector<double> &y,
		   std::vector<double> &q, std::vector<double> &L);
    double stepsize_select(const std::vector<double> &y);
    bool passes_system_test(std::vector<double> &y);

    void set_gas_data_ptr(Gas_data &Q)
    { Q_ = &Q; }
    
    int get_n_reactions()
    { return reaction_.size(); }
    
    Reaction* get_reaction( int ir )
    { return reaction_[ir]; }

    int get_directional_rates( std::vector<double> &w_f, std::vector<double> &w_b );

private:
    // A list of Reactions making up the reaction scheme
    std::vector<Reaction*> reaction_;
    double err_tol_;
    // A ragged array of integers. The list corresponds
    // to the number of species.  Each entry in the list
    // is a list itself, which denotes which reactions
    // a given species participates in.
    std::vector<std::vector<int> > participation_;
    Gas_data *Q_;
    std::vector<double> q_, L_, ydot_;
    std::vector<double> massf_, M_;
};

#endif
