// Author: Daniel F Potter
// Date: 07-Dec-2009

#ifndef MT_NONEQUILIBRIUM_HH
#define MT_NONEQUILIBRIUM_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "reaction-rate-coeff.hh"
#include "../models/gas_data.hh"
#include "generalised-Arrhenius.hh"
#include "../models/species-energy-modes.hh"

class MarroneTreanor_dissociation : public Generalised_Arrhenius {
public:
    MarroneTreanor_dissociation(lua_State *L, Gas_model &g, double T_upper, double T_lower);
    MarroneTreanor_dissociation(double A, double n, double E_a, double T_upper, double T_lower,
				double U, std::string v_name);
    ~MarroneTreanor_dissociation();

private:
    double U_;
    std::vector<Species_energy_mode*> vib_modes_;
    int iTv_;
    
private:
    int s_eval(const Gas_data &Q);
};

Reaction_rate_coefficient* create_MarroneTreanor_dissociation_coefficient(lua_State *L, Gas_model &g, double T_upper, double T_lower);

#endif
