// Author: Rowan J. Gollan
// Date: 21-Oct-2008
// Place: Poquoson, Virginia, USA
//

#ifndef EQUILIBRIUM_CONSTANT_HH
#define EQUILIBRIUM_CONSTANT_HH

#include <map>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../models/gas-model.hh"


class Equilibrium_constant {
public:
    Equilibrium_constant( int iT );
    virtual ~Equilibrium_constant() {};
    
    double eval(const Gas_data &Q)
    { return s_eval(Q); }
    
    int get_iT()
    { return iT_; }

protected:
    int iT_;
    
protected:
    virtual double s_eval(const Gas_data &Q) = 0;
    
};

Equilibrium_constant* create_Equilibrium_constant(lua_State *L,
						  std::map<int, int> &nu,
						  Gas_model &g);

#endif
