// \file fluid_thermo.cxx
//
// Thermodynamic properties for a pure fluid.
// Since our interest in mainly in the gas phase for the following code, 
// the phase transitions have not been dealt with explicitly.
//
// Adapted from the text
// W C Reynolds
// Thermodynamic properties in SI: graphs, tables and 
// computational equations for 40 substances.
//
// PJ, Apr 2008

#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "fluid_thermo.hh"

// using namespace std;
using std::cout;
using std::cerr;
using std::endl;

// string representation of state
std::string State::str(int precision) const
{
    std::ostringstream ost;
    ost << std::setprecision(precision);
    ost << "p=" << p/1.0e6 << "MPa T=" << T << "K h=" << h/1.0e3 
	<< "kJ/kg s=" << s/1.0e3 << "kJ/kg.K rho=" << rho << "kg/m**3 u="
	<< u/1.0e3 << "kJ/kg";
    return ost.str();
}

PureFluid::PureFluid( std::string name )
    : name(name)
{}

PureFluid::~PureFluid() {}

double PureFluid::p_rhoT( double rho, double T )
{
    return 0.0;
}

double PureFluid::dpdT_rhoT( double rho, double T )
{
    // Should put in a finite-difference calculation here.
    return 0.0;
}

double PureFluid::p_sat( double T )
{
    return 0.0;
}

double PureFluid::dpsatdT( double T )
{
    return 0.0;
}

double PureFluid::sat_liquid_density( double T )
{
    return 0.0;
}

double PureFluid::cv0( double T )
{
    return 0.0;
}

int PureFluid::us_rhoT( double rho, double T, double &u, double &h, double &s )
{
    // Integrate the equation of state to get internal energy and entropy.
    double sum_u = u0;
    double sum_s = s0 - R*log(rho);
    // First part of integral path (zero density, ideal gas).
    int n1 = 100;
    double dT = (T - T0) / n1;
    double Tmid, cv;
    for ( int i = 0; i < n1; ++i ) {
	Tmid = T0 + (0.5+i)*dT;
	cv = cv0(Tmid);
	sum_u += cv * dT;
	sum_s += cv/Tmid * dT;
    }
    // Second part of integral path (constant T).
    int n2 = 100;
    double drho = rho / n2;
    double rhomid, p, dpdT;
    for ( int i = 0; i < n2; ++i ) {
	rhomid = (0.5+i)*drho;
	p = p_rhoT(rhomid, T);
	dpdT = dpdT_rhoT(rhomid, T);
	sum_u += drho * (p - T*dpdT) / (rhomid*rhomid);
	sum_s += drho * (rhomid*R - dpdT) / (rhomid*rhomid);
    }
    p = p_rhoT(rho, T);
    u = sum_u;
    h = sum_u + p/rho;
    s = sum_s;
    return 0;
}

int PureFluid::eval_state( State *st, std::string option,
			   int initial_guess_supplied )
{
    // Utility function to get all thermodynamic properties given two
    // as indicated by the option string.
    // option  values given
    // ---------------------------------
    // "DT"    density and temperature
    // "PT"    pressure and temperature
    // "ST"    entropy and temperature
    // "PH"    pressure and enthalpy
    // "PS"    pressure and entropy
    // "DU"    density and internal energy
    //
    // I had grand plans for this code to be very neat but that's all gone to pot...
    // Also, had to use a pointer to State st because of SWIG and python.
    //
    double rho_given, T_given, p_given, u_given, h_given, s_given;
    double delta_rho, delta_T;
    // Tolerances for terminating the Newton iterations are selected
    // depending on the epected range of the particular variable.
    double p_tolerance = 1.0;
    double s_tolerance = 0.001;
    double h_tolerance = 1.0;
    double u_tolerance = 1.0;
    double f1, f2, df1drho, df2drho, df1dT, df2dT, drho, dT, denom;
    int count;
    // The following are temporary variables for us in evaluating finite differences.
    // _01 represents a step in the rho coordinate
    // _10 represents a step in the T coordinate
    double rho_00, rho_01, T_00, T_10;
    double p_00, p_01, p_10, u_00, u_01, u_10;
    double h_00, h_01, h_10, s_00, s_01, s_10;

    if ( option == "dt" || option == "DT" ) {
	// density and temperature are given, drop through to final evaluation
    } else if ( option == "pt" || option == "PT" ) {
	p_given = st->p;
	T_given = st->T;
	// We need to obtain the correct density.
	if ( initial_guess_supplied ) {
	    rho_00 = st->rho;
	} else {
	    rho_00 = 1.0; // somewhere to start when we know nothing
	}
	T_00 = T_given;
	p_00 = p_rhoT(rho_00, T_00);
	f1 = p_00 - p_given;
	count = 0;
	while ( fabs(f1) > p_tolerance && count < 20 ) {
	    drho = 0.001 * rho_00 + 0.0001;
	    rho_01 = rho_00 + drho;
	    p_01 = p_rhoT(rho_01, T_00);
	    df1drho = (p_01 - p_00) / drho;
	    delta_rho = -f1 / df1drho;
	    rho_00 += delta_rho;
	    T_00 = T_given;
	    p_00 = p_rhoT(rho_00, T_00);
	    f1 = p_00 - p_given;
	    ++count;
	}
	st->rho = rho_00;
	st->T = T_given;
	if ( fabs(f1) > p_tolerance ) {
	    cout << "state evaluation: iteration to find p failed, count=" 
		 << count << endl;
	}
    } else if ( option == "st" || option == "ST" ) {
	s_given = st->s;
	T_given = st->T;
	// We need to obtain the correct density.
	if ( initial_guess_supplied ) {
	    rho_00 = st->rho;
	} else {
	    rho_00 = 1.0; // somewhere to start when we know nothing
	}
	T_00 = T_given;
	us_rhoT(rho_00, T_00, u_00, h_00, s_00);
	f1 = s_00 - s_given;
	count = 0;
	while ( fabs(f1) > s_tolerance && count < 20 ) {
	    drho = 0.001 * rho_00 + 0.0001;
	    rho_01 = rho_00 + drho;
	    us_rhoT(rho_01, T_00, u_01, h_01, s_01);
	    df1drho = (s_01 - s_00) / drho;
	    delta_rho = -f1 / df1drho;
	    rho_00 += delta_rho;
	    T_00 = T_given;
	    us_rhoT(rho_00, T_00, u_00, h_00, s_00);
	    f1 = s_00 - s_given;
	    ++count;
	}
	st->rho = rho_00;
	st->T = T_given;
	if ( fabs(f1) > s_tolerance ) {
	    cout << "state evaluation: iteration to find s failed, count=" 
		 << count << endl;
	}
    } else if ( option == "ph" || option == "PH" ) {
	p_given = st->p;
	h_given = st->h;
	// We need to obtain the correct density.
	if ( initial_guess_supplied ) {
	    rho_00 = st->rho;
	    T_00 = st->T;
	} else {
	    rho_00 = 1.0; // somewhere to start when we know nothing
	    T_00 = 400.0;
	}
	p_00 = p_rhoT(rho_00, T_00);
	f1 = p_00 - p_given;
	us_rhoT(rho_00, T_00, u_00, h_00, s_00);
	f2 = h_00 - h_given;
	count = 0;
	while ( (fabs(f1) > p_tolerance || fabs(f2) > h_tolerance) && count < 20 ) {
	    drho = 0.001 * rho_00 + 0.0001;
	    rho_01 = rho_00 + drho;
	    p_01 = p_rhoT(rho_01, T_00);
	    df1drho = (p_01 - p_00) / drho;
	    us_rhoT(rho_01, T_00, u_01, h_01, s_01);
	    df2drho = (h_01 - h_00) / drho;

	    dT = 0.001 * T_00 + 0.01;
	    T_10 = T_00 + dT;
	    p_10 = p_rhoT(rho_00, T_10);
	    df1dT = (p_10 - p_00) / dT;
	    us_rhoT(rho_00, T_10, u_10, h_10, s_10);
	    df2dT = (h_10 - h_00) / dT;

	    denom = (df1drho*df2dT - df1dT*df2drho);
	    delta_rho = (df1dT*f2 - df2dT*f1) / denom;
	    rho_00 += delta_rho;
	    delta_T = (df2drho*f1 - df1drho*f2) / denom;
	    T_00 += delta_T;
	    p_00 = p_rhoT(rho_00, T_00);
	    f1 = p_00 - p_given;
	    us_rhoT(rho_00, T_00, u_00, h_00, s_00);
	    f2 = h_00 - h_given;
	    ++count;
	}
	st->rho = rho_00;
	st->T = T_00;
	if ( fabs(f1) > p_tolerance || fabs(f2) > h_tolerance ) {
	    cout << "state evaluation: iteration to find p,h failed, count=" 
		 << count << endl;
	}
    } else if ( option == "ps" || option == "PS" ) {
	p_given = st->p;
	s_given = st->s;
	// We need to obtain the correct density.
	if ( initial_guess_supplied ) {
	    rho_00 = st->rho;
	    T_00 = st->T;
	} else {
	    rho_00 = 1.0; // somewhere to start when we know nothing
	    T_00 = 400.0;
	}
	p_00 = p_rhoT(rho_00, T_00);
	f1 = p_00 - p_given;
	us_rhoT(rho_00, T_00, u_00, h_00, s_00);
	f2 = s_00 - s_given;
	count = 0;
	while ( (fabs(f1) > p_tolerance || fabs(f2) > s_tolerance) && count < 20 ) {
	    drho = 0.001 * rho_00 + 0.0001;
	    rho_01 = rho_00 + drho;
	    p_01 = p_rhoT(rho_01, T_00);
	    df1drho = (p_01 - p_00) / drho;
	    us_rhoT(rho_01, T_00, u_01, h_01, s_01);
	    df2drho = (s_01 - s_00) / drho;

	    dT = 0.001 * T_00 + 0.01;
	    T_10 = T_00 + dT;
	    p_10 = p_rhoT(rho_00, T_10);
	    df1dT = (p_10 - p_00) / dT;
	    us_rhoT(rho_00, T_10, u_10, h_10, s_10);
	    df2dT = (s_10 - s_00) / dT;

	    denom = (df1drho*df2dT - df1dT*df2drho);
	    delta_rho = (df1dT*f2 - df2dT*f1) / denom;
	    delta_rho = (fabs(delta_rho) > 5.0) ? delta_rho/3.0 : delta_rho;
	    rho_00 += delta_rho;
	    delta_T = (df2drho*f1 - df1drho*f2) / denom;
	    delta_T = (fabs(delta_T) > 500.0) ? delta_T/3.0 : delta_T;
	    T_00 += delta_T;
	    p_00 = p_rhoT(rho_00, T_00);
	    f1 = p_00 - p_given;
	    us_rhoT(rho_00, T_00, u_00, h_00, s_00);
	    f2 = s_00 - s_given;
	    ++count;
	    // cout << "count=" << count << " rho=" << rho_00 << " T=" 
	    //      << T_00 << " p=" << p_00 << " s=" << s_00 << endl;
	}
	st->rho = rho_00;
	st->T = T_00;
	if ( fabs(f1) > p_tolerance || fabs(f2) > s_tolerance ) {
	    cout << "state evaluation: iteration to find p,s failed, count=" 
		 << count << endl;
	}

    } else if ( option == "du" || option == "DU" ) {
	rho_given = st->rho;
	u_given = st->u;
	// We need to obtain the correct density.
	if ( initial_guess_supplied ) {
	    T_00 = st->T;
	} else {
	    T_00 = 400.0; // somewhere to start when we know nothing
	}
	rho_00 = rho_given;
	us_rhoT(rho_00, T_00, u_00, h_00, s_00);
	f1 = u_00 - u_given;
	count = 0;
	while ( fabs(f1) > u_tolerance && count < 20 ) {
	    dT = 0.001 * T_00 + 0.01;
	    T_10 = T_00 + dT;
	    us_rhoT(rho_00, T_10, u_10, h_10, s_10);
	    df1dT = (u_10 - u_00) / dT;
	    delta_T = -f1 / df1dT;
	    T_00 += delta_T;
	    rho_00 = rho_given;
	    us_rhoT(rho_00, T_00, u_00, h_00, s_00);
	    f1 = u_00 - u_given;
	    ++count;
	}
	st->T = T_00;
	st->rho = rho_given;
	if ( fabs(f1) > u_tolerance ) {
	    cout << "state evaluation: iteration to find u failed, count=" 
		 << count << endl;
	}
    } else {
	cerr << "Error: eval_state(): unknown option " << option << endl;
	return 1;
    }
    // At this point we should have density and temperature correct.
    // Fill in the other properties consistently.
    st->p = p_rhoT(st->rho, st->T);
    us_rhoT(st->rho, st->T, st->u, st->h, st->s);
    return 0;
}

//-----------------------------------------------------------------------
// Helper functions that may be reused for several fluids.
// These essentially implement WC Reynolds' equations 
// from the appendices of his text.

double eqnP3( double rho, double T, 
	      double R, double gama, std::vector<double> &A )
{
    // Equation P-3 for pressure, given density and temperature.
    // As transcribed from Reynold's text and manipulated a little.
    double Tinv = 1.0/T;
    double p = rho * R * T;
    double rhos = rho * rho;
    p += rhos * (A[1]*T + A[2] + Tinv*(A[3] + Tinv*(A[4] + Tinv*A[5])));
    rhos *= rho;
    p += rhos * (T*A[6] + A[7] + Tinv*A[8]);
    rhos *= rho;
    p += rhos * (T*A[9] + A[10]);
    rhos *= rho;
    p += rhos * (T*A[11] + A[12]);
    rhos *= rho;
    p += rhos * A[13];
    // last two terms
    rhos = rho * rho * rho;
    double term1 = rhos * Tinv*Tinv*(A[14] + Tinv*(A[15] + Tinv*A[16]));
    rhos *= rho * rho;
    double term2 = rhos * Tinv*Tinv*(A[17] + Tinv*(A[18] + Tinv*A[19]));
    p += (term1 + term2) * exp(-gama*rho*rho);
    return p;
}

double ddTeqnP3( double rho, double T, 
		 double R, double gama, std::vector<double> &A )
{
    // Derivative (with respect to T) of equation P-3 for pressure, 
    // given density and temperature.
    // This is needed for the integrals for internal energy, and entropy.
    double Tinv = 1.0/T;
    double dpdT = rho * R;
    double rhos = rho * rho;
    dpdT += rhos * (A[1] - Tinv*Tinv*(A[3] + Tinv*(2.0*A[4] + Tinv*3.0*A[5])));
    rhos *= rho;
    dpdT += rhos * (A[6] - Tinv*Tinv*A[8]);
    rhos *= rho;
    dpdT += rhos * A[9];
    rhos *= rho;
    dpdT += rhos * A[11];
    // last two terms
    rhos = rho * rho * rho;
    double term1 = (-1.0)*rhos*Tinv*Tinv*Tinv*(2.0*A[14]+Tinv*(3.0*A[15]+Tinv*4.0*A[16]));
    rhos *= rho * rho;
    double term2 = (-1.0)*rhos*Tinv*Tinv*Tinv*(2.0*A[17]+Tinv*(3.0*A[18]+Tinv*4.0*A[19]));
    dpdT += (term1 + term2) * exp(-gama*rho*rho);
    return dpdT;
}

double eqnS2( double T, double Tc, double Tp, double Pc_sat, 
	      std::vector<double> &F )
{
    double TTp = T/Tp - 1.0;
    double sum = F[1]+TTp*(F[2]+TTp*(F[3]+TTp*(F[4]+TTp*(F[5]+TTp*(F[6]+TTp*(F[7]+TTp*F[8]))))));
    double rhs = (Tc/T - 1.0) * sum;
    double p_sat = exp(rhs) * Pc_sat;
    return p_sat;
}

double ddTeqnS2( double T, double Tc, double Tp, double Pc_sat, 
		 std::vector<double> &F )
{
    double p_sat = eqnS2(T, Tc, Tp, Pc_sat, F);
    double dT = 1.0;
    if ( T > (Tc-2.0) ) dT = -1.0;
    double p2 = eqnS2(T+dT, Tc, Tp, Pc_sat, F);
    // simple one-sided difference
    double D1 = (p2 - p_sat)/dT;
    // do it again with a smaller step size
    dT *= 0.5;
    p2 = eqnS2(T+dT, Tc, Tp, Pc_sat, F);
    double D2 = (p2 - p_sat)/dT;
    // extrapolate the difference to limit of dT=0
    return 2.0*D2 - D1;
}

double eqnD2( double T, double Tc, std::vector<double> &D )
{
    // Equation D2 for saturated liquid density, given temperature.
    double X = 1.0 - T/Tc;
    double rho_f = D[1] + D[2]*pow(X,1.0/3) + D[3]*pow(X,2.0/3) +
	D[4]*X + D[5]*pow(X,4.0/3) + D[6]*pow(X,5.0/3);
    return rho_f;
}

double eqnC6( double T, std::vector<double> &G )
{
    // Equation C6 for specific heat, const volume, ideal gas.
    double c_v_0 = G[1]/T + G[2] + T*(G[3] + T*(G[4]+ T*(G[5] + T*G[6])));
    return c_v_0;
}

//------------------------------------------------------------------
// Definitions for specific fluids...

CarbonDioxide::CarbonDioxide()
    : PureFluid("carbon-dioxide")
{
    M = 44.01;       // kg/kmol
    Tc = 304.21;     // degree K
    Pc = 7.3834e6;   // Pa
    Pc_sat = 7.38350e6; // Pa
    rhoc = 464.00;   // kg/m**3
    T0 = 216.54;     // degree K
    R = 188.918;     // J/kg.K
    A.resize(20);
    A[0] =  0.0;
    A[1] =  2.2488558e-1;
    A[2] = -1.3717965e+2;
    A[3] = -1.4430214e+4;
    A[4] = -2.9630491e+6;
    A[5] = -2.0606039e+8;
    A[6] =  4.5554393e-5;
    A[7] =  7.7042840e-2;
    A[8] =  4.0602371e+1;
    A[9] =  4.0029509e-7;
    A[10]= -3.9436077e-4;
    A[11]=  1.2115286e-10;
    A[12]=  1.0783386e-7;
    A[13]=  4.3962336e-11;
    A[14]= -3.6505545e+4;
    A[15]=  1.9490511e+7;
    A[16]= -2.9186718e+9;
    A[17]=  2.4358627e-2;
    A[18]= -3.7546530e+1;
    A[19]=  1.1898141e+4;
    gama = 5.0e-6;
    D.resize(7);
    D[0] =  0.0;
    D[1] =  4.6400009e+2;
    D[2] =  6.7938129e+2;
    D[3] =  1.4776836e+3;
    D[4] = -3.1267676e+3;
    D[5] =  3.6397656e+3;
    D[6] = -1.3437098e+3;
    F.resize(9);
    F[0] =  0.0;
    F[1] = -6.5412610;
    F[2] = -2.7914636e-1;
    F[3] = -3.4716202;
    F[4] = -3.4989637;
    F[5] = -1.9770948e+1;
    F[6] =  1.3922839e+2;
    F[7] = -2.7670389e+2;
    F[8] = -7.0510251e+3;
    Tp = 250.0;
    G.resize(7);
    G[0] =  0.0;
    G[1] =  8.726361e+3;
    G[2] =  1.840040e+2;
    G[3] =  1.914025;
    G[4] = -1.667825e-3;
    G[5] =  7.305950e-7;
    G[6] = -1.255290e-10;
    u0 = 3.2174105e+5;
    s0 = 2.1396056e+3;
}

CarbonDioxide::~CarbonDioxide()
{
    A.resize(0);
    D.resize(0);
    F.resize(0);
    G.resize(0);
}

double CarbonDioxide::p_rhoT( double rho, double T )
{
    return eqnP3( rho, T, R, gama, A );
}

double CarbonDioxide::dpdT_rhoT( double rho, double T )
{
    return ddTeqnP3( rho, T, R, gama, A );
}

double CarbonDioxide::p_sat( double T )
{
    return eqnS2( T, Tc, Tp, Pc_sat, F );
}

double CarbonDioxide::dpsatdT( double T )
{
    return ddTeqnS2( T, Tc, Tp, Pc_sat, F );
}

double CarbonDioxide::sat_liquid_density( double T )
{
    return eqnD2( T, Tc, D );
}

double CarbonDioxide::cv0( double T )
{
    return eqnC6( T, G );
}
