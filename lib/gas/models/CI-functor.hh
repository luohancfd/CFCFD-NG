// Author: Daniel F Potter
// Version: 12-Oct-2009
//          Initial coding.
// Version: 13-Apr-2010
//          Changed name from CI_functor to GuptaYos_CI_functor.

#ifndef CI_FUNCTOR_HH
#define CI_FUNCTOR_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/functor.hh"

#include "chemical-species.hh"

/******** Gupta-Yos curve fit ********/

class GuptaYos_CI_functor : public Univariate_functor {
public:
    GuptaYos_CI_functor( std::string section, double T_low, double T_high, lua_State *L );
    GuptaYos_CI_functor( const GuptaYos_CI_functor &c );
    ~GuptaYos_CI_functor();
    
    GuptaYos_CI_functor* clone() const;
    
    double operator()(double T);
    
private:
    double T_low_;
    double T_high_;
    std::vector<double> a_;
};

/******** Bruno curve fits ********/

// Heavy particle interaction (neutral/neutral or neutral/ion)
class Bruno_HPI_CI_functor : public Univariate_functor {
public:
    Bruno_HPI_CI_functor( int l, int s, int Z_l, int Z_s, lua_State * L );
    Bruno_HPI_CI_functor( const Bruno_HPI_CI_functor &c );
    ~Bruno_HPI_CI_functor();
    
    Bruno_HPI_CI_functor* clone() const;
    
    double operator()(double T);
    
private:
    double eps0_;
    double sigma2_;
    std::vector<double> a_;
};

// Charged particle interaction
class Bruno_CPI_CI_functor : public Bivariate_functor {
public:
    Bruno_CPI_CI_functor( int l, int m, lua_State * L );
    Bruno_CPI_CI_functor( const Bruno_CPI_CI_functor &c );
    ~Bruno_CPI_CI_functor();
    
    Bruno_CPI_CI_functor* clone() const;
    
    double operator()(double T, double N_e);
    
private:
    double eps0_;
    double constA_;
    double constB_;
};

/******** Wright curve fits ********/

// Charged particle interaction - shielded Coloumb potential
class Stallcop_SCP_CI_functor : public Bivariate_functor {
public:
    Stallcop_SCP_CI_functor( int l, int m, int Z_l, int Z_m );
    Stallcop_SCP_CI_functor( const Stallcop_SCP_CI_functor &c );
    ~Stallcop_SCP_CI_functor();
    
    Stallcop_SCP_CI_functor* clone() const;
    
    double operator()(double T, double N_e);
    
private:
    double C_;
    double c_;
    double D_;
};

/********** Neufeld curve fits ***********/

// Empirical equations to calculate CI's for the Leonard Jones (12-6) potential
class Neufeld_CI_functor : public Univariate_functor {
public:
    Neufeld_CI_functor( int l, int s, Chemical_species * I, Chemical_species * J, lua_State * L );
    Neufeld_CI_functor( const Neufeld_CI_functor &n );
    ~Neufeld_CI_functor();
    
    Neufeld_CI_functor* clone() const;
    
    double operator()(double T);
    
private:
    double eps0_K_;
    double sigma_;
    std::vector<double> c_;

};

#endif
