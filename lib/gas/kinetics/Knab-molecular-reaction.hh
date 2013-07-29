// Author: Daniel F Potter
// Date: 07-Dec-2009

#ifndef KNAB_NONEQUILIBRIUM_HH
#define KNAB_NONEQUILIBRIUM_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "reaction-rate-coeff.hh"
#include "../models/gas_data.hh"
#include "generalised-Arrhenius.hh"
#include "../models/species-energy-modes.hh"

class Knab_molecular_reaction : public Generalised_Arrhenius {
public:
    Knab_molecular_reaction(lua_State *L, Gas_model &g);
    Knab_molecular_reaction(double A, double n, double E_a, double U, double alpha, std::string v_name);
    ~Knab_molecular_reaction();

private:
    double U_;
    double alpha_;
    double alpha_A_;
    std::vector<Species_energy_mode*> vib_modes_;
    int iTv_;
    
private:
    int s_eval(const Gas_data &Q);
};

Reaction_rate_coefficient* create_Knab_molecular_reaction_coefficient(lua_State *L, Gas_model &g);

#endif
