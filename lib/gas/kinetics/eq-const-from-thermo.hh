// Author: Rowan J. Gollan
// Date: 21-Oct-2008
// Place: Poquoson, Virginia, USA
//

#ifndef EQ_CONST_FROM_THERMO
#define EQ_CONST_FROM_THERMO

#include <map>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "equilibrium-constant.hh"
#include "../models/gas-model.hh"


class Eq_const_from_thermo : public Equilibrium_constant {
public:
    Eq_const_from_thermo(std::map<int, int> &nu, Gas_model &g, int iT);
    ~Eq_const_from_thermo();
private:
    double s_eval(const Gas_data &Q);
    std::map<int, int> nu_;
    Gas_model &g_;
    Gas_data * Q_;
};

Equilibrium_constant* create_Eq_const_from_thermo(lua_State *L,
						  std::map<int, int> &nu,
						  Gas_model &g);

#endif
