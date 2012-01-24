// Author: Rowan J. Gollan
// Date: 04-Nov-2008
// Place: Hampton, Virginia, USA

#ifndef CEA_CP_FUNCTOR_HH
#define CEA_CP_FUNCTOR_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/functor.hh"

class CEA_Cp_functor : public Univariate_functor {
public:
    CEA_Cp_functor(lua_State *L, double R);
    CEA_Cp_functor(const CEA_Cp_functor &c);
    ~CEA_Cp_functor();
    CEA_Cp_functor* clone() const;

    double operator()(double T);

private:
    double R_;
    double T_low_;
    double T_high_;
    std::vector<double> a_;
};

#endif
