// Author: Rowan J. Gollan
// Date: 12-Sep-2008

#ifndef REACTION_UPDATE_HH
#define REACTION_UPDATE_HH

#include <string>
#include <vector>

#include "../models/gas_data.hh"
#include "../models/gas-model.hh"

class Reaction_update {
public:
    Reaction_update() {}
    virtual ~Reaction_update() {}
    
    int update_state(Gas_data &Q, double dt, double &dt_suggest, Gas_model *gm=0)
    { return s_update_state(Q, dt, dt_suggest, gm); }

    int rate_of_change(Gas_data &Q, std::vector<double> &dcdt)
    { return s_rate_of_change(Q, dcdt); }

    double update_state_py(Gas_data &Q, double dt, double dt_suggest=-1.0, Gas_model *gm=0)
    { s_update_state(Q, dt, dt_suggest, gm); return dt_suggest; }

    std::vector<double> rate_of_change_py(Gas_data &Q)
    { 
	std::vector<double> dcdt(Q.massf.size(), 0.0);
	s_rate_of_change(Q, dcdt); 
	return dcdt; 
    }
    
    int eval_chemistry_energy_coupling_source_terms( Gas_data &Q, std::vector<double> &dedt )
    { return s_eval_chemistry_energy_coupling_source_terms( Q, dedt ); }

    int get_directional_rates( std::vector<double> &w_f, std::vector<double> &w_b )
    { return s_get_directional_rates( w_f, w_b ); }

protected:
    virtual int s_update_state(Gas_data &Q, double dt, double &dt_suggest, Gas_model *gm) = 0;
    virtual int s_rate_of_change(Gas_data &Q, std::vector<double> &dcdt) = 0;
    virtual int s_eval_chemistry_energy_coupling_source_terms( Gas_data &Q, std::vector<double> &dedt ) = 0;
    virtual int s_get_directional_rates( std::vector<double> &w_f, std::vector<double> &w_b ) = 0;

};

Reaction_update* create_Reaction_update(std::string cfile, Gas_model &g);

#endif
