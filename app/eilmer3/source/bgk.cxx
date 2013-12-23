/// \file bgk.cxx
/// \ingroup eilmer3
/// \brief General functions to do with the solution of non-equilibrium flows. 
///
/// \author DB
/// \version 14-Sep-12 initial coding

#include <math.h>

#include "bgk.hh"


/// \brief Reduced form of the model Shakhov equation
Vector3 Shakhov(double rho, double U, double V, double T, double qx, double qy, double R, double Pr, double u, double v)
{
    // pre-compute powers and common factors
    double RT = R*T;
    double RT3 = RT*RT*RT;

    double uU = u - U;
    double vV = v - V;
    
    double uU2 = uU*uU;
    double vV2 = vV*vV;
    double UV = uU2 + vV2;

    double EXP =rho/(2.0*M_PI*exp(UV/(2.0*RT)));
    double Q = (Pr - 1)*(qx*uU + qy*vV);
    double SUB = 5*RT3*rho;

    Vector3 gh;

    gh.x = (1/RT)*EXP*(1 - (Q*(UV-4*RT))/SUB);
    gh.y = EXP*(1 - (Q*(UV-2*RT))/SUB);

    return gh;
}
