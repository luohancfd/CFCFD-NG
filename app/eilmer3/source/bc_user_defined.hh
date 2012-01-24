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
    int nsp, nmodes;
    lua_State *L;
public:
    UserDefinedBC( Block &bdp, int which_boundary, 
		   const std::string filename="udf.lua",
		   bool is_wall=false, bool use_udf_flux=false );
    UserDefinedBC( const UserDefinedBC &bc );
    virtual ~UserDefinedBC();
    virtual int apply_inviscid( double t ); // copies flow data to ghost cells
    virtual int apply_viscous( double t ); // sets wall T to user-defined value
private:
    int eval_flux_udf( double t, int i, int j, int k, FV_Interface *IFace );
    int eval_inviscid_udf( double t, int i, int j, int k, FV_Interface *IFace );
    CFlowCondition *unpack_flow_table( void );
    int eval_viscous_udf( double t, int i, int j, int k, FV_Interface *IFace, FV_Cell *cell );
    void handle_lua_error(lua_State *L, const char *fmt, ...);
};

//--------------------------------------------------------------------
// Functions that are registered with the Lua interpreter have to be
// outside the class definitions.

int luafn_sample_flow(lua_State *L);
int luafn_locate_cell(lua_State *L);
int luafn_compute_diffusion_coefficient(lua_State *L);

//--------------------------------------------------------------------

class AdjacentPlusUDFBC : public UserDefinedBC {
public:
    AdjacentPlusUDFBC( Block &bdp, int which_boundary, int other_block, 
		       int other_face, int neighbour_orientation=0,
		       const std::string filename="udf.lua",
		       bool is_wall=false, bool use_udf_flux=false );
    AdjacentPlusUDFBC( const AdjacentPlusUDFBC &bc );
    virtual ~AdjacentPlusUDFBC();
};
