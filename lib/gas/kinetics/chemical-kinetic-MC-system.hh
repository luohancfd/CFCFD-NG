// Author: Rowan J. Gollan
// Date: 12-Sep-2008

#ifndef CHEMICAL_KINETIC_MC_SYSTEM_HH
#define CHEMICAL_KINETIC_MC_SYSTEM_HH

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
#include "species-pieces.hh"
#include "chemistry-energy-coupling.hh"

class Chemical_kinetic_MC_system : public OdeSystem {
public:
    Chemical_kinetic_MC_system(lua_State *L, Gas_model &g, int nreac, double error_tol,
			       double T_upper, double T_lower);
    Chemical_kinetic_MC_system(std::string cfile, Gas_model &g, int nreac);

    ~Chemical_kinetic_MC_system();
    
    int eval(const std::vector<double> &y, std::vector<double> &ydot);
    int eval_new_concentrations(const std::vector<double> &y, std::vector<double> &c);
    int eval_species_rates(const std::vector<double> &y, std::vector<double> &cdot);
    double stepsize_select(const std::vector<double> &y);
    bool passes_system_test(std::vector<double> &y);
    void print_reaction_rates();
    void print_species_rates();
    void print_limiting_species_and_reaction();

    void set_gas_data_ptr_and_initial_concs(Gas_data &Q, std::vector<double> &c);
    
    int get_n_reactions()
    { return reaction_.size(); }
    
    Reaction* get_reaction( int ir )
    { return reaction_[ir]; }
    
    int initialise_chemistry_energy_coupling( Gas_data &Q, std::vector<double> &c_old );
    
    int apply_chemistry_energy_coupling( Gas_data &Q, std::vector<double> &delta_c,
    				         std::vector<double> &c_new );
    
    int eval_chemistry_energy_coupling_source_terms( Gas_data &Q, const std::vector<double> &y, std::vector<double> &dedt );

    size_t cecs_size( void ) { return cecs_.size(); }

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
    std::vector<double> ydot_;
    std::vector<double> w_;
    std::vector<double> c_;
    std::vector<double> cinit_, massf_, M_;
    // A list of 'Species_pieces' to manage concentration changes
    std::vector<Species_pieces*> spec_;
    // A list of 'Energy_coupling' (one for each reaction as needed) 
    // to manage energy changes due to reactions
    std::vector<Chemistry_energy_coupling*> cecs_;
};

#endif
