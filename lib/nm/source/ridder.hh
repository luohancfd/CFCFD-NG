// ridder.hh

#ifndef RIDDER_HH
#define RIDDER_HH

#include <functional>

double solve(std::function<double (double)>f, double x1, double x2, 
	     double tol=1.0e-9);
int bracket(std::function<double (double)>f, double &x1, double &x2,
	    int max_try=50, double factor=1.6);

#endif
