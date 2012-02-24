// Author: Daniel F Potter
// Version: 12-Oct-2009
//          Initial coding.
// Version: 10-Dec-2009
//          Split into transport and diffusion flavours
// Version: 13-Apr-2010
//          Complete overhaul - this class now just computes Omega_ij and D_ij

#ifndef COLLISION_INTEGRAL_HH
#define COLLISION_INTEGRAL_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/segmented-functor.hh"

#include "gas_data.hh"
#include "chemical-species.hh"

class Collision_integral {
public:
    Collision_integral() {}
    Collision_integral( int iT, std::string type );
    virtual ~Collision_integral() {}
    
    double eval_Pi_Omega_11( Gas_data &Q )
    { return s_eval_Pi_Omega_11(Q); }
    
    double eval_Pi_Omega_22( Gas_data &Q )
    { return s_eval_Pi_Omega_22(Q); }
    
    // Delete this entry eventually
    double eval_D( Gas_data &Q )
    { return s_eval_D(Q); }
    
    std::string get_type()
    { return type_; }

protected:
    virtual double s_eval_Pi_Omega_11( Gas_data &Q ) = 0;
    virtual double s_eval_Pi_Omega_22( Gas_data &Q ) = 0;
    virtual double s_eval_D( Gas_data &Q ) = 0;
    
    int iT_;
    std::string type_;
};

class No_CI_model : public Collision_integral {
public:
    No_CI_model() {}
    ~No_CI_model() {}
    
private:
    double s_eval_Pi_Omega_11( Gas_data &Q )
    { return 0.0; }
    
    double s_eval_Pi_Omega_22( Gas_data &Q )
    { return 0.0; }
    
    double s_eval_D( Gas_data &Q )
    { return 0.0; }
};

class GuptaYos_CI_model : public Collision_integral {
public:
    GuptaYos_CI_model( int iT, int Z_i, int Z_j, lua_State * L );
    ~GuptaYos_CI_model();
    
private:
    double s_eval_Pi_Omega_11( Gas_data &Q );
    double s_eval_Pi_Omega_22( Gas_data &Q );	
    double s_eval_D( Gas_data &Q );

    int Z_i_;
    int Z_j_;
    Segmented_functor * Pi_Omega_11_;
    Segmented_functor * Pi_Omega_22_;
    Segmented_functor * D_;
};

class Stallcop_CI_model : public Collision_integral {
public:
    Stallcop_CI_model( int iT, int Z_i, int Z_j, lua_State * L );
    ~Stallcop_CI_model();
    
private:
    double s_eval_Pi_Omega_11( Gas_data &Q );
    double s_eval_Pi_Omega_22( Gas_data &Q );	
    double s_eval_D( Gas_data &Q );
    
    int ie_;
    
    Bivariate_functor * Pi_Omega_11_;
    Bivariate_functor * Pi_Omega_22_;
};

class Bruno_CI_model : public Collision_integral {
public:
    Bruno_CI_model( int iT, int Z_i, int Z_j, lua_State * L );
    ~Bruno_CI_model();
    
private:
    double s_eval_Pi_Omega_11( Gas_data &Q );
    double s_eval_Pi_Omega_22( Gas_data &Q );	
    double s_eval_D( Gas_data &Q );
    
    Univariate_functor * Pi_Omega_11_;
    Univariate_functor * Pi_Omega_22_;
};

class Neufeld_CI_model : public Collision_integral {
public:
    Neufeld_CI_model( int iT, Chemical_species * I, Chemical_species * J, lua_State * L );
    ~Neufeld_CI_model();
    
private:
    double s_eval_Pi_Omega_11( Gas_data &Q );
    double s_eval_Pi_Omega_22( Gas_data &Q );	
    double s_eval_D( Gas_data &Q );
    
    Univariate_functor * Pi_Omega_11_;
    Univariate_functor * Pi_Omega_22_;
};

Collision_integral * new_CI_from_file( std::string fname );

#endif
