// Author: Rowan J. Gollan
// Date: 10-Oct-2008

#ifndef GENERALISED_ARRHENIUS_HH
#define GENERALISED_ARRHENIUS_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "reaction-rate-coeff.hh"
#include "../models/gas_data.hh"

class Generalised_Arrhenius : public Reaction_rate_coefficient {
public:
    Generalised_Arrhenius(lua_State *L, Gas_model &g);
    Generalised_Arrhenius(double A, double n, double E_a);
    ~Generalised_Arrhenius();
    
    double get_A() { return A_; }
    double get_n() { return n_; }
    double get_E_a() { return E_a_; }
    
private:
    int s_eval(const Gas_data &Q);
    int s_eval_from_T(const double T);
    double A_;
    double n_;
    double E_a_;
};

Reaction_rate_coefficient* create_Generalised_Arrhenius_coefficient(lua_State *L, Gas_model &g);

#endif
