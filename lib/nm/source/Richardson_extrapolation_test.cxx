// Author: Rowan J. Gollan
// Date: 08-Jul-2008

#include <cstdio>
#include <cmath>

#include "functor.hh"
#include "Richardson_extrapolation.hh"

class Test_functor : public Univariate_functor {
    double operator()(double x)
    {
	return pow(x, 2.0)*cos(x);
    }
};

int main()
{
    int i, j, status;
    double x, h, deriv;
    Test_functor f;

    h = 0.1;
    x = 1.0;

    deriv = R_extrap_deriv(f, x, h, status);

    printf("Testing Richardson extrapolation to estimate derivative.\n");
    printf("Numerical answer: %.8f\n", deriv);
    printf("Exact answer    : %.8f\n", 2.0*cos(1.0) - sin(1.0));
    
    return 0;
}

	
	
