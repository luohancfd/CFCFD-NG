// Author: Daniel F Potter
// Version: 13-Apr-2010
//          Initial coding.

#ifndef BINARY_INTERACTION_HH
#define BINARY_INTERACTION_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "collision-integral.hh"

class Binary_interaction {
public:
    Binary_interaction( int isp, int jsp, lua_State * L );
    ~Binary_interaction();
    
    void store_Delta_1( Gas_data &Q )
    { Delta_1_ = s_eval_Delta_1(Q); }
    
    void store_Delta_2( Gas_data &Q )
    { Delta_2_ = s_eval_Delta_2(Q); }
    
    double get_Delta_1()
    { return Delta_1_; }
    
    void set_Delta_1( double Delta_1 )
    { Delta_1_ = Delta_1; }
    
    double get_Delta_2()
    { return Delta_2_; }
    
    void set_Delta_2( double Delta_2 )
    { Delta_2_ = Delta_2; }
    
    double eval_D( Gas_data &Q )
    { return s_eval_D(Q); }
    
    void store_D( Gas_data &Q )
    { D_val_ = s_eval_D(Q); }
    
    double get_D()
    { return D_val_; }
    
    void set_D( double D )
    { D_val_ = D; }
    
    Collision_integral * get_CI_model_ptr()
    { return CI_model_; }
    
public:
    int isp_;
    int jsp_;
    
private:
    double s_eval_Delta_1( Gas_data &Q );
    double s_eval_Delta_2( Gas_data &Q );
    double s_eval_D( Gas_data &Q );
    
    double m_i_;
    double m_j_;
    int iT_;
    Collision_integral * CI_model_;
    
    double Delta_1_;
    double Delta_2_;
    double D_val_;
};

#endif
