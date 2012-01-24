/** \file secant.hh
 *  \ingroup nm
 *  \brief Declarations for the secant function solver.
 *  \author PJ
 *  \version 01-Jan-2006
 */

#ifndef SECANT_HH
#define SECANT_HH

#include <string>
#include <iostream>
#include <vector>
using namespace std;

double secant_solve( double (*f)(double), 
		     double x0, double x1,
		     int &result_flag,
		     double tolerance=1.0e-11, 
		     int max_iterations=100 );

#endif
