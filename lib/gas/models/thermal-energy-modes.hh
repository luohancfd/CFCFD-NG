// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#ifndef THERMAL_ENERGY_MODE_HH
#define THERMAL_ENERGY_MODE_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "species-energy-modes.hh"
#include "chemical-species.hh"

class Thermal_energy_mode {
public:
    Thermal_energy_mode( std::string name, std::vector<Chemical_species*> &species, lua_State *L );
    virtual ~Thermal_energy_mode();

    int no_components()
    { return component_names_.size(); };

    std::string component_name(int ic)
    { return component_names_[ic]; }
    
    double decode_conserved_energy(Gas_data &Q, double rhoe)
    { return s_decode_conserved_energy(Q,rhoe); }
    
    double encode_conserved_energy(const Gas_data &Q)
    { return s_encode_conserved_energy(Q); }
    
    // NOTE: following functions give values per kg of this mode, i.e.
    //       NOT yet weighted by how much fraction this mode occupies
    //       of total
    //       08-Jan-2012 : Change by R. Gollan
    //                     Earlier implementation gave J/kg-mixture
    //                     (so the value was already weighted)
    double eval_energy(const Gas_data &Q)
    { return s_eval_energy(Q); }
    
    double eval_temperature(Gas_data &Q)
    { return s_eval_temperature(Q); }
    
    double eval_Cv( Gas_data &Q )
    { return s_eval_dedT( Q ); }
    
    void test_derivatives( Gas_data &Q )
    { return s_test_derivatives(Q); }

    double mode_massf(const Gas_data &Q);
    
    std::string get_name()
    { return s_get_name(); }

protected:
    int iT_;
    
    std::string name_;
    std::vector<int> sp_idx_; // List of species indices for compute mass fraction in mode
    std::vector<std::string> component_names_;
    std::vector<Species_energy_mode*> components_;
    
    std::string s_get_name()
    { return name_; }

    virtual double s_decode_conserved_energy(Gas_data &Q, double rhoe);
    virtual double s_encode_conserved_energy(const Gas_data &Q);
    virtual double s_eval_dedT(Gas_data &Q);
    virtual double s_eval_energy( const Gas_data &Q);
    virtual double s_eval_temperature(Gas_data &Q) = 0;
    virtual void s_test_derivatives(Gas_data &Q) = 0;
};

class Constant_Cv_energy_mode : public Thermal_energy_mode {
public:
    Constant_Cv_energy_mode( std::string name, std::vector<Chemical_species*> &species, lua_State *L );
    ~Constant_Cv_energy_mode() {}
    
private:
    double s_eval_temperature( Gas_data &Q );
    void s_test_derivatives(Gas_data &Q);
};

class Variable_Cv_energy_mode : public Thermal_energy_mode {
public:
    Variable_Cv_energy_mode( std::string name, std::vector<Chemical_species*> &species, lua_State *L );
    ~Variable_Cv_energy_mode() {}
    
private:
    double T_min_, T_max_;
    int max_iterations_;
    double convergence_tolerance_;
    
    bool check_T_range( double T );
    void impose_T_limits( double &T );
    bool check_convergence( double e, double e_given );
    double f_zero( double e, double e_given );
    double dfdT( Gas_data &Q ) { return s_eval_dedT( Q ); }
    
    double s_eval_temperature(Gas_data &Q);
    double s_eval_temperature_bisection(Gas_data &Q, double toll);
    void s_test_derivatives(Gas_data &Q);
};

#endif
