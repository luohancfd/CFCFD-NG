// Author: Daniel F. Potter
// Date: 02-Dec-2009

#ifndef ENERGY_EXCHANGE_UPDATE_HH
#define ENERGY_EXCHANGE_UPDATE_HH

#include <string>
#include <vector>

#include "../models/gas_data.hh"
#include "../models/gas-model.hh"

class Energy_exchange_update {
public:
    Energy_exchange_update() {}
    virtual ~Energy_exchange_update() {}
    
    int update_state(Gas_data &Q, double dt, double &dt_suggest, Gas_model *gm)
    { return s_update_state(Q, dt, dt_suggest, gm); }

    int rate_of_change(Gas_data &Q, std::vector<double> &dedt)
    { return s_rate_of_change(Q, dedt); }

    double update_state_py(Gas_data &Q, double dt, double dt_suggest=-1.0, Gas_model *gm=0)
    { s_update_state(Q, dt, dt_suggest, gm); return dt_suggest; }

    std::vector<double> rate_of_change_py(Gas_data &Q)
    { 
	std::vector<double> dedt(Q.massf.size(), 0.0);
	s_rate_of_change(Q, dedt); 
	return dedt; 
    }

protected:
    virtual int s_update_state(Gas_data &Q, double dt, double &dt_suggest, Gas_model *gm) = 0;
    virtual int s_rate_of_change(Gas_data &Q, std::vector<double> &dedt) = 0;

};

Energy_exchange_update* create_Energy_exchange_update(std::string cfile, Gas_model &g);

#endif
