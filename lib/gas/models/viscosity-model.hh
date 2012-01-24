// Author: Rowan J. Gollan
// Version: 24-May-2008
//            Initial coding.
//

#ifndef VISCOSITY_MODEL_HH
#define VISCOSITY_MODEL_HH

#include "gas_data.hh"

class Viscosity_model {
public:
    Viscosity_model() {}
    virtual ~Viscosity_model() {}

    double eval_viscosity(const Gas_data &Q)
    { return s_eval_viscosity(Q); }

private:
    virtual double s_eval_viscosity(const Gas_data &Q) = 0;
};

#endif
