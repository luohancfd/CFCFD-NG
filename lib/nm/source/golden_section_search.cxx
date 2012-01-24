// Author: Rowan J. Gollan
// Version: 06-Jun-2008
//            Initial coding

#include <cmath>

#include "../../util/source/useful.h"

double golden_section_search(double (*f)(double),
			     double a, double c,
			     int &result_flag,
			     double tolerance,
			     int max_iterations)
{
    double g = (3.0 - sqrt(5.0))/2.0;
    double b1 = a + g*(c - a);
    double b2 = a + (1.0 - g)*(c - a);
    double f1 = (*f)(b1);
    double f2 = (*f)(b2);

    if( (c - a) < tolerance*(a + c) ) {
	result_flag = SUCCESS;
	return (a + c)/2.0;
    }

    for( int i = 0; i < max_iterations; ++i ) {
	if( f1 <= f2) {
	    c = b2;
	    b2 = b1;
	    f2 = f1;
	    b1 = a + g*(c - a);
	    f1 = (*f)(b1);
	}
	else {
	    a = b1;
	    b1 = b2;
	    f1 = f2;
	    b2 = a + (1 - g)*(c - a);
	    f2 = (*f)(b2);
	}
	
	// Test if we are close enough
	if( fabs(c - a) <= tolerance ) {
	    result_flag = SUCCESS;
	    return (a + c)/2.0;
	}
    }

    // If we make it this far we've exceed max_iterations
    result_flag = ITERATION_ERROR;
    return (a + c)/2.0;
		
}
