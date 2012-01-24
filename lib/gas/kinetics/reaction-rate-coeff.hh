// Author: Rowan J. Gollan
// Version: 10-Oct-2008
//          (Port of rr_coeffs.hh/.cxx to new implementation)

#ifndef REACTION_RATE_COEFF_HH
#define REACTION_RATE_COEFF_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../models/gas-model.hh"

class Reaction_rate_coefficient {
public:
    Reaction_rate_coefficient() {};
    virtual ~Reaction_rate_coefficient() {};
    
    double k()
    { return k_; }

    int eval(const Gas_data &Q)
    { return s_eval(Q); }
    
    int eval_from_T(const double T)
    { return s_eval_from_T(T); }

protected:
    double k_;
    virtual int s_eval(const Gas_data &Q) = 0;
    virtual int s_eval_from_T(const double T);
};

#ifndef SWIG
Reaction_rate_coefficient* create_Reaction_rate_coefficient(lua_State *L, Gas_model &g);
#endif
Reaction_rate_coefficient* create_Reaction_rate_coefficient(std::string cfile, std::string rate, Gas_model &g);

#endif
