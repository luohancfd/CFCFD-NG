#ifndef FLOW_STATE_HH
#define FLOW_STATE_HH

#include <string>
#include <vector>

#include "../../../lib/gas/models/gas_data.hh"

class Flow_state {
public:
    Flow_state();
    Flow_state(Gas_data &Q1, double u1, double Q_rad1=0.0);
    Flow_state(Flow_state &fs); 
    ~Flow_state();

    void set_flow_state(Gas_data &Q1, double u1, double Q_rad1=0.0);

    std::string str(bool with_Q_rad);
    std::string str(bool with_Q_rad, std::string specied_output_type, const std::vector<double> &M);
    Gas_data * Q;
    double u;
    double Q_rad;
};

#endif
