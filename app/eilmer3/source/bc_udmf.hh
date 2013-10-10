#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../../lib/gas/models/gas-model.hh"
#include "block.hh"
#include "bc.hh"

class UserDefinedMassFluxBC : public BoundaryCondition {
public:
    UserDefinedMassFluxBC(Block *bdp, int which_boundary, 
			  const std::string fname="udf.lua");
    UserDefinedMassFluxBC(const UserDefinedMassFluxBC &bc);
    UserDefinedMassFluxBC();
    UserDefinedMassFluxBC& operator=(const UserDefinedMassFluxBC &bc);
    virtual ~UserDefinedMassFluxBC();
    void print_info(std::string lead_in);
    // default apply_convective() is just to reflect normal velocity
    int apply_viscous(double t); // sets wall T to user-defined value

private:
    std::string filename_;
    lua_State *L;
    Gas_model *gmodel;
    size_t nsp, nmodes;
    int start_interpreter();
    int eval_user_fn(double t, size_t i, size_t j, size_t k, FV_Interface *IFace);
    int handle_lua_error(lua_State *L, const char *fmt, ...);
};

