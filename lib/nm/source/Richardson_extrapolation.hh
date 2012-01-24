// Author: Rowan J. Gollan
// Date: 08-Jul-2008

#ifndef RICHARDSON_EXTRAPOLATION_HH
#define RICHARDSON_EXTRAPOLATION_HH

#include "functor.hh"

double R_extrap_deriv(Univariate_functor &f, double x, double h, int &status,
		      double tolerance=1.0e-6, int max_steps=5);

#endif
