// Author: Rowan J. Gollan
// Date: 08-Jul-2008

#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>

#include "../../util/source/useful.h"
#include "Richardson_extrapolation.hh"

using namespace std;

double R_extrap_deriv(Univariate_functor &f, double x, double h, int &status,
		      double tolerance, int max_steps)
{
    int i, j;
    vector<vector<double> > d;
    
    // Maximum rows in table..
    d.resize(max_steps);
    
    // First row.
    d[0].resize(1);
    double d0 = (f(x + h) - f(x - h))/(2*h);

    d[0][0] = d0;

    // Begin building table
    for ( i = 1; i < max_steps; ++i ) {
	h = h / 2.0;
	d[i].resize(i+1);
	d[i][0] = (f(x + h) - f(x - h))/(2*h);
	for ( j = 1; j <= i; ++j ) {
	    // Apply extrapolation to fill out row
	    d[i][j] = d[i][j-1] + (d[i][j-1] - d[i-1][j-1])/(pow(2.0, 2.0*j) - 1.0);
	}

	// Check for convergence
	if ( fabs(d[i][i] - d[i][i-1]) <= tolerance ) {
	    status = SUCCESS;
	    return d[i][i];
	}
    }

    // If we get this far, we've failed.
    status = FAILURE;
    return d[i-1][i-1];
}


