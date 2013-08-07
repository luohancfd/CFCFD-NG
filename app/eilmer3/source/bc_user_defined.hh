// bc_user_defined.hh


extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}
#include "bc.hh"

class UserDefinedBC : public BoundaryCondition {
public:
    std::string filename;
private:
    CFlowCondition *fdata1;
    CFlowCondition *fdata2;
    Gas_model *gmodel;
    size_t nsp, nmodes;
    lua_State *L;
public:
    UserDefinedBC(Block *bdp, int which_boundary, 
		  const std::string filename="udf.lua",
		  bool is_wall=false);
    UserDefinedBC(const UserDefinedBC &bc);
    UserDefinedBC();
    UserDefinedBC & operator=(const UserDefinedBC &bc);
    virtual ~UserDefinedBC();
    virtual int apply_convective(double t); // copies flow data to ghost cells
    virtual int apply_viscous(double t); // sets wall T to user-defined value
private:
    int start_interpreter();
    int eval_conv_flux_udf(double t, size_t i, size_t j, size_t k, FV_Interface *IFace);
    int eval_ghost_cell_udf(double t, size_t i, size_t j, size_t k, FV_Interface *IFace);
    CFlowCondition *unpack_flow_table(void);
    int eval_iface_udf(double t, size_t i, size_t j, size_t k, FV_Interface *IFace, const FV_Cell *cell);
    int eval_visc_flux_udf(double t, size_t i, size_t j, size_t k, FV_Interface *IFace);
    void handle_lua_error(lua_State *L, const char *fmt, ...);
};

//--------------------------------------------------------------------

class AdjacentPlusUDFBC : public UserDefinedBC {
public:
    AdjacentPlusUDFBC(Block *bdp, int which_boundary, int other_block, 
		      int other_face, int neighbour_orientation=0,
		      const std::string filename="udf.lua",
		      bool is_wall=false);
    AdjacentPlusUDFBC(const AdjacentPlusUDFBC &bc);
    virtual ~AdjacentPlusUDFBC();
};
