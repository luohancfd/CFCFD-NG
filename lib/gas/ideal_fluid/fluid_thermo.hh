// \file fluid_thermo.hh
//
// Thermodynamic properties for a pure fluid.
//
// Adapted from the text
// W C Reynolds
// Thermodynamic properties in SI: graphs, tables and 
// computational equations for 40 substances.
//
// PJ, Apr 2008

#ifndef FLUID_THERMO_HH
#define FLUID_THERMO_HH

#include <string>
#include <vector>

// The following struct will be useful for passing information 
// back to the Python interpreter (when we use it).
class State {
public:
    double p;   // pressure, Pa
    double T;   // temperature, degrees K
    double rho; // density kg/m**3
    double u;   // internal energy, J/kg
    double h;   // enthalpy, J/kg
    double s;   // entropy, J/kg/K
    std::string str(int precision=4) const; // string representation
};

class PureFluid {
public:
    std::string name;

public:
    double Tc, Pc, rhoc; // critical conditions
    double Pc_sat;       // may be slightly different in the 
    double M;            // molecular mass
    double T0;           // triple-point temperature
    double R;            // gas constant
    double u0;           // reference internal energy
    double s0;           // reference entropy

public:
    PureFluid( std::string name );
    virtual ~PureFluid();
    virtual double p_rhoT( double rho, double T ); 
    virtual double dpdT_rhoT( double rho, double T );
    virtual double p_sat( double T );
    virtual double dpsatdT( double T );
    virtual double sat_liquid_density( double T );
    virtual double cv0( double T );
    virtual int us_rhoT( double rho, double T, double &u, double &h, double &s );
    virtual int eval_state( State *st, std::string option="dt",
			    int initial_guess_supplied=0 );
};

// Helper functions that may be reused for several fluids.
double eqnP3( double rho, double T, 
	      double R, double gama, std::vector<double> &A );
double ddTeqnP3( double rho, double T, 
		 double R, double gama, std::vector<double> &A );
double eqnS2( double T, double Tc, double Tp, double Pc_sat, 
	      std::vector<double> &F );
double ddTeqnS2( double T, double Tc, double Tp, double Pc_sat, 
		 std::vector<double> &F );
double eqnD2( double T, double Tc, std::vector<double> &D );
double eqnC6( double T, std::vector<double> &G );

class CarbonDioxide : public PureFluid {
public:
    CarbonDioxide();
    virtual ~CarbonDioxide();
    virtual double p_rhoT( double rho, double T );
    virtual double dpdT_rhoT( double rho, double T );
    virtual double p_sat( double T );
    virtual double dpsatdT( double T );
    virtual double sat_liquid_density( double T );
    virtual double cv0( double T );
protected:
    // Constants for Reynold's curve fits.
    std::vector<double> A, D, F, G;
    double Tp, gama;
};

#endif
