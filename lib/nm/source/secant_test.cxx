/** \file secant_test.cxx
 *  \ingroup nm
 *  \brief Test the secant-method function solver. 
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
#include "secant.hh"

// Example from Gerald and Wheatley, p. 45
// solution: x=1.732051
double f( double x ) { return ((x + 1.0) * x - 3.0) * x - 3.0; }

int main()
{
    double x0, x1, x2;
    int result_flag;
    cout << "Begin demonstration of secant_solve()..." << endl;
    x0 = 1.0; x1 = 2.0; 
    x2 = secant_solve(f, x0, x1, result_flag);
    cout << "Solution at x=" << x2 << endl;
    cout << "Check: f(" << x2 << ")=" << f(x2) << endl;
    cout << "Done." << endl;
    return 0;
}
