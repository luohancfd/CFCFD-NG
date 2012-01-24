/** \file secant.cxx
 *  \ingroup nm
 *  \brief Implementation of the secant-method function solver. 
 *  \author PJ
 *  \version 01-Jan-2006 initial coding
 *
 * This module is a C++ replacement for the combined Python zero_solver code 
 * that was built by Rowan.
 */
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#include "../../util/source/useful.h"
#include "secant.hh"

double secant_solve( double (*f)(double), 
		     double x0, double x1,
		     int &result_flag,
		     double tolerance,
		     int max_iterations )
{
    double x2, f0, f1, df, temp;
    int i;
    // Evaluate the initial guesses.
    f0 = (*f)(x0);
    f1 = (*f)(x1);

    if ( fabs(f0) < fabs(f1) ) {
	// Ensure that x1 is the closer guess.
        temp = x0; x0 = x1; x1 = temp;
        temp = f0; f0 = f1; f1 = temp;
    }
    if ( fabs(f1) < tolerance ) {
	result_flag = SUCCESS;
	return x1;
    }
    // Now, repeatedly improve the guess.
    for ( i = 0; i < max_iterations; ++i) {
	df = f0 - f1;
	if ( df == 0.0 ) {
	    x2 = 0.5 * (x0 + x1);  // Maybe this will help.
	} else {
	    x2 = x1 - f1 * (x0 - x1) / df; // Secant-method update.
	}
	x0 = x1; f0 = f1; 
	x1 = x2; f1 = (*f)(x2); // current guess and its evaluation
        if ( fabs(f1) < tolerance ) break;
    }
    if ( i >= (max_iterations-1) ) {
	cout << "The secant method did not converge after " 
	     << i+1 << " iterations." << endl;
	result_flag = ITERATION_ERROR;
	return x1;
    }
    result_flag = SUCCESS;
    return x1;
}
