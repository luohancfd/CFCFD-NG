// Author: Daniel F. Potter
// Date: 21-Sep-2009
// Place: Brisbane, Queendland, AUST

#ifndef NONEQ_THERMAL_BEHAVIOUR_HH
#define NONEQ_THERMAL_BEHAVIOUR_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/segmented-functor.hh"
#include "gas_data.hh"
#include "thermal-behaviour-model.hh"
#include "thermal-energy-modes.hh"
#include "chemical-species.hh"
#include "chemical-equilibrium-system.hh"

class Noneq_thermal_behaviour : public Thermal_behaviour_model {
public:
    Noneq_thermal_behaviour(lua_State *L);
    ~Noneq_thermal_behaviour();
    
    int get_number_of_modes()
    { return (int) modes_.size(); }

    int mode_no_components(int imode)
    { return modes_[imode]->no_components(); }
    
    std::string mode_name(int imode)
    { return modes_[imode]->get_name(); }

    std::string mode_component_name(int imode, int ic)
    { return modes_[imode]->component_name(ic); }

    int eval_equilibrium_composition( double T, double rho, std::vector<double> &massf )
    { return ces_->solve_system(T,rho,massf); }
    
    int test_chemical_equilibrium_system( double T, double p, std::vector<double> &molef )
    { return ces_->test_system(T,p,molef); }
    
    Partial_equilibrium_reaction * get_partial_equilibrium_reaction_pointer( size_t index )
    { return ces_->get_partial_equilibrium_reaction_pointer( index ); }
    
    void test_derivatives_for_mode( int itm, Gas_data &Q )
    { return modes_[itm]->test_derivatives(Q); }

private:
    double min_massf_;
    std::vector<Thermal_energy_mode*> modes_;
    std::vector<Chemical_species*> species_;
    // For equilibrium composition determination
    Chemical_equilibrium_system * ces_;

    int s_decode_conserved_energy(Gas_data &Q, const std::vector<double> &rhoe);
    int s_encode_conserved_energy(const Gas_data &Q, std::vector<double> &rhoe);
    double s_dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status);
    double s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status);
    int s_eval_energy(Gas_data &Q, Equation_of_state *EOS_);
    int s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_);
    double s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_modal_enthalpy_isp( const Gas_data &Q, Equation_of_state *EOS_, int isp, int itm );
    double s_eval_modal_Cv(Gas_data &Q, Equation_of_state *EOS_, int itm );
    double s_eval_modal_massf(const Gas_data &Q, int itm);
};

Noneq_thermal_behaviour * new_ntb_from_file( std::string inFile );

#endif
