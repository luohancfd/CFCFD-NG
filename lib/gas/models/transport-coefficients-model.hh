// Author: Rowan J. Gollan
// Date: 09-Jul-2008

#ifndef TRANSPORT_COEFFICIENTS_MODEL_HH
#define TRANSPORT_COEFFICIENTS_MODEL_HH

#include "gas_data.hh"

class Transport_coefficients_model {
public:
    Transport_coefficients_model() {}
    virtual ~Transport_coefficients_model() {}

    int eval_transport_coefficients(Gas_data &Q)
    { return s_eval_transport_coefficients(Q); }

private:
    virtual int s_eval_transport_coefficients(Gas_data &Q) = 0;
};

#endif
