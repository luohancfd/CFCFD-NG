// Author: Rowan J. Gollan
// Date: 04-Nov-2008
// Place: Hampton, Virginia, USA

#ifndef CEA_H_FUNCTOR_HH
#define CEA_H_FUNCTOR_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/functor.hh"

class CEA_h_functor : public Univariate_functor {
public:
    CEA_h_functor(lua_State *L, double Cp_low, double Cp_high, double R);
    CEA_h_functor(const CEA_h_functor &h);
    ~CEA_h_functor();
    CEA_h_functor* clone() const;

    double operator()(double T);

private:
    double eval(double T);
    double R_;
    double T_low_;
    double T_high_;
    double Cp_low_;
    double Cp_high_;
    std::vector<double> a_;
};

#endif
