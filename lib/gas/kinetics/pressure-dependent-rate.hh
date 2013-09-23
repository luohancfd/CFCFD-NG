// Author: Rowan J. Gollan
// Date: 16-Apr-2009
// Place: NIA, Hampton, Virginia, USA
//
// This a port from Brendan O'Flaherty's implementation
// found in gas_models2.
//

#ifndef PRESSURE_DEPENDENT_RATE_HH
#define PRESSURE_DEPENDENT_RATE_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "reaction-rate-coeff.hh"
#include "generalised-Arrhenius.hh"
#include "../models/gas_data.hh"

double compute_third_body_value(const Gas_data &Q, std::map<int, double> efficiencies, std::vector<double> M);

class Pressure_dependent : public Reaction_rate_coefficient {
public:
    Pressure_dependent(lua_State *L, Gas_model &g, double T_upper, double T_lower);
    ~Pressure_dependent();
    
    double get_third_body_value(const Gas_data &Q) { return compute_third_body_concentration(Q); }

private:
    int s_eval(const Gas_data &Q);
    double compute_third_body_concentration(const Gas_data &Q) { return compute_third_body_value(Q, efficiencies_, M_); }

    std::vector<double> M_;

    std::map<int, double> efficiencies_;

    Generalised_Arrhenius *k_inf_;
    Generalised_Arrhenius *k_0_;

    // Some (possible) Troe model parameters
    bool Troe_model_;
    bool T2_supplied_;
    double a_;
    double T1_;
    double T2_;
    double T3_;

};

Reaction_rate_coefficient* create_pressure_dependent_coefficient(lua_State *L, Gas_model &g,
								 double T_upper, double T_lower);

#endif
