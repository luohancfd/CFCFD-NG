#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../../lib/gas/models/gas-model.hh"
#include "block.hh"
#include "bc.hh"

class UserDefinedEnergyFluxBC : public BoundaryCondition {
public:
    UserDefinedEnergyFluxBC(Block *bdp, int which_boundary, 
			  const std::string fname="udf.lua");
    UserDefinedEnergyFluxBC(const UserDefinedEnergyFluxBC &bc);
    UserDefinedEnergyFluxBC();
    UserDefinedEnergyFluxBC& operator=(const UserDefinedEnergyFluxBC &bc);
    virtual ~UserDefinedEnergyFluxBC();
    void print_info(std::string lead_in);
    // default apply_convective() is just to reflect normal velocity
    int apply_viscous(double t); 

private:
    std::string filename_;
    lua_State *L;
    Gas_model *gmodel;
    size_t nsp, nmodes;
    int start_interpreter();
    int eval_user_fn(double t, size_t i, size_t j, size_t k, FV_Interface *IFace);
    int handle_lua_error(lua_State *L, const char *fmt, ...);
};

