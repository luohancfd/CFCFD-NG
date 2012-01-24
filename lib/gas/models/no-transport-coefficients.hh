// Author: Daniel. F Potter
// Date: 24-Sep-2009

#ifndef NO_TRANSPORT_COEFFICIENTS_HH
#define NO_TRANSPORT_COEFFICIENTS_HH

#include "gas_data.hh"
#include "transport-coefficients-model.hh"

class No_transport_coefficients : public Transport_coefficients_model {
public:
    No_transport_coefficients() {}
    ~No_transport_coefficients() {}

private:
    int s_eval_transport_coefficients(Gas_data &Q);
};

#endif
