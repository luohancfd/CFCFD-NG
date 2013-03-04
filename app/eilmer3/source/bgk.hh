/** \file bgk.hh
 * \ingroup eilmer3
 * \brief Header file containing general functions to do with the solution of non-equilibrium flows.
 * 
 * \author DB
 * \version 14-Sep-12 initial coding
 */


#ifndef BGK_HH
#define BGK_HH

#include "../../../lib/geometry2/source/geom.hh"

const int PI = 3.141592653589793;

Vector3 Shakhov(double rho, double U, double V, double T, 
		double qx, double q, double R, double Pr,
		double u, double v);

#endif
