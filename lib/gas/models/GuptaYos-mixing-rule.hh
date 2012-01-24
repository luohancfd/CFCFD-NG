// Author: Daniel F Potter
// Date: 13-Oct-2009

#ifndef GUPTAYOS_MIXING_RULE_HH
#define GUPTAYOS_MIXING_RULE_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "transport-coefficients-model.hh"
#include "binary-interaction.hh"
#include "chemical-species.hh"

class GuptaYos_mixing_rule : public Transport_coefficients_model {
public:
    GuptaYos_mixing_rule(lua_State *L);
    ~GuptaYos_mixing_rule();
    
    Binary_interaction * get_binary_interaction_ptr( int isp, int jsp );
    Collision_integral * get_collision_integral_ptr( int isp, int jsp );

private:
    std::vector<Binary_interaction*> unique_BIs_;
    std::vector< std::vector<Binary_interaction*> > BI_table_;

    int nsp_, nsp_se_, e_index_;
    double ignore_mole_fraction_;
    std::vector<double> m_;
    std::vector<int> Z_;
    std::vector<double> x_;
    std::vector<Chemical_species*> species_;

    int s_eval_transport_coefficients(Gas_data &Q);
};

GuptaYos_mixing_rule * create_GuptaYos_mixing_rule_from_file( std::string lua_file );

#endif
