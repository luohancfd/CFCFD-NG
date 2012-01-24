// Author: Rowan J. Gollan
// Version: 24-May-2008
//            Initial coding.
//

#ifndef THERMAL_CONDUCTIVITY_MODEL_HH
#define THERMAL_CONDUCTIVITY_MODEL_HH

#include "gas_data.hh"

class Thermal_conductivity_model {
public:
    Thermal_conductivity_model() {}
    virtual ~Thermal_conductivity_model() {}

    double eval_thermal_conductivity(const Gas_data &Q)
    { return s_eval_thermal_conductivity(Q); }

private:
    virtual double s_eval_thermal_conductivity(const Gas_data &Q) = 0;
};

#endif
