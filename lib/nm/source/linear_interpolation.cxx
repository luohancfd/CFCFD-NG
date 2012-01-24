// Authors: RJG and BTO
// Date: 22-Nov-2009

#include <iostream>

#include "../../util/source/useful.h"
#include "../../util/source/dbc_assert.hh"
#include "linear_interpolation.hh"

using namespace std;

int 
linear_eval(double xval,
	    double &yval,
	    const vector<double> x,
	    const vector<double> y)
{
    int i0, i1;
    double wA, wB;

    // ASSERT(x.size() == y.size());
    
    // trivial cases
    if (x.size() == 0) {
	return FAILURE;
    }
    if (x.size() == 1) {
	yval = y[0];
	return SUCCESS;
    }

    // Implementation of a 1-D lookup-table
    // with linear interpolation

    // Search for appropriate index (lower-bound)
    i0 = 0;

    for ( ; i0 < (int)x.size(); ++i0) {
	if ( x[i0] >= xval ) {
	    i0--;
	    break;
	}
    }
    
    // Establish weight values

    // 1. Handle edges of table
    if( i0 == -1 ) {
	i0 = 0;
	i1 = i0 + 1;
	wA = 1.0;
	wB = 0.0;
    }
    else if( i0 == int(x.size())) {
	i0 = int(x.size()-2);
	i1 = i0 + 1;
	wA = 0.0;
	wB = 1.0;
    }
    // 2. and all other points (interior)
    else {
	i1 = i0 + 1;
	wA = (x[i1] - xval)/(x[i1] - x[i0]);
	wB = (xval - x[i0])/(x[i1] - x[i0]);
    }

    // Perform the interpolation once the
    // weights are established.
    yval = wA*y[i0] + wB*y[i1];
    
    return SUCCESS;
}

double 
linear_eval_py(double xval,
	       std::vector<double> x,
	       std::vector<double> y) 
{
    double yval; 
    linear_eval(xval, yval, x, y); 
    return yval; 
}
