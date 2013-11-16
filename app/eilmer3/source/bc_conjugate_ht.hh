// bc_conjugate_ht.hh
#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "bc.hh"


class ConjugateHeatTransferBC : public BoundaryCondition {
public:
    ConjugateHeatTransferBC(Block *bdp, int which_boundary)
    ConjugateHeatTransferBC(const ConjugateHeatTransferBC &bc);
    ConjugateHeatTransferBC();
    ConjugateHeatTransferBC& operator=(const ConjugateHeatTransferBC &bc);
    virtual ~ConjugateHeatTransferBC();
    void print_info(std::string lead_in);
    // default apply_convective() -- just want to reflect normal velocity
    int apply_viscous(double t); // uses T from wall model
};
