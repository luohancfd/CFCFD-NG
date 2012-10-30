// Author: Daniel F. Potter
// Date: 30-Oct-2012
// Place: Gšttingen, Germany
//

#ifndef EQ_CONST_FROM_PARTITION_FUNCTIONS
#define EQ_CONST_FROM_PARTITION_FUNCTIONS

#include <map>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "equilibrium-constant.hh"
#include "../models/gas-model.hh"


class Eq_const_from_partition_functions : public Equilibrium_constant {
public:
    Eq_const_from_partition_functions(std::map<int, int> &nu, Gas_model &g, int iT);
    ~Eq_const_from_partition_functions();
private:
    double s_eval(const Gas_data &Q);
    std::map<int, int> nu_;
    Gas_model &g_;
};

Equilibrium_constant* create_Eq_const_from_partition_functions(lua_State *L,
						  std::map<int, int> &nu,
						  Gas_model &g);

#endif
