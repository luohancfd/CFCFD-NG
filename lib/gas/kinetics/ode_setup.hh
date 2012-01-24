// Author: Rowan J. Gollan
// Date: 17-Oct-2008
// Place: Hampton, Virginia, USA

#ifndef ODE_SETUP_HH
#define ODE_SETUP_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/ode_solver.hh"

OdeSolver* create_ode_solver(lua_State *L, int isp, std::string name="anonymous");

#endif
