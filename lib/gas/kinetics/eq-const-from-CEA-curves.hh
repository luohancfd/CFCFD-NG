// Author: Daniel F. Potter
// Date: 20-Apr-2010
// Place: Dutton Park, Brisbane, QLD
//

#ifndef EQ_CONST_FROM_CEA_CURVES
#define EQ_CONST_FROM_CEA_CURVES

#include <map>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "equilibrium-constant.hh"
#include "../models/gas-model.hh"


class Eq_const_from_CEA_curves : public Equilibrium_constant {
public:
    Eq_const_from_CEA_curves(std::map<int, int> &nu, Gas_model &g, int iT);
    ~Eq_const_from_CEA_curves();
private:
    double s_eval(const Gas_data &Q);
    std::map<int, int> nu_;
    Gas_model &g_;
};

Equilibrium_constant* create_Eq_const_from_CEA_curves(lua_State *L, std::map<int, int> &nu, Gas_model &g);

#endif
