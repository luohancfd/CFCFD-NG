// Author: Rowan J. Gollan
// Date: 09-Jul-2008

#ifndef TRANSPORT_COEFFICIENTS_MODEL_HH
#define TRANSPORT_COEFFICIENTS_MODEL_HH

#include "gas_data.hh"
#include "gas-model.hh"
#include "thermal-behaviour-model.hh"

class Transport_coefficients_model {
public:
    Transport_coefficients_model() {}
    virtual ~Transport_coefficients_model() {}

    int eval_transport_coefficients(Gas_data &Q, Gas_model *gmodel=0)
    { return s_eval_transport_coefficients(Q, gmodel); }

private:
    virtual int s_eval_transport_coefficients(Gas_data &Q, Gas_model *gmodel=0) = 0;
};

#endif
