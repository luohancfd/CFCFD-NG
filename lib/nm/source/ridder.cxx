// ridder.cxx
// Solve a nonlinear equation f(x)=0 using the method of Ridder.
// Peter J. 
// mech3750 demo code 12-Mar-2014
// added bracketing 19-Mar-2014

#include <stdexcept>
#include <math.h>
#include "ridder.hh"

double solve(std::function<double (double)>f, double x1, double x2, double tol) 
{
    // Locate a root of f(x) by subdividing the original range,
    // assuming a linear model of the underlying transformed function.
    //
    // Input:
    //    f: user-supplied function f(x)
    //    x1: first end of range
    //    x2: other end of range
    // Returns:
    //    x, a point near the root.
    //
    double x3, f3;
    double x4 = x1; // So that g++ doesn't warn on maybe unitialized.
    double f4, eps;

    double f1 = f(x1); 
    double f2 = f(x2);
    if ( fabs(f1) == 0.0 ) return x1;
    if ( fabs(f2) == 0.0 ) return x2;
    if ( x1 == x2 ) {
	throw std::runtime_error("Bad initial range given to bracket.");
    }
    if ( f1 * f2 > 0.0 ) {
	throw std::runtime_error("Range does not clearly bracket a root.");
    }
    while ( fabs(x2 - x1) > tol ) {
	x3 = 0.5*(x1+x2);
	f3 = f(x3);
	if ( f3 == 0.0 ) return x3;
	eps = (f3 + copysign(sqrt(f3*f3-f1*f2),f2))/f2;
	x4 = x3 - f3*eps*(x1-x3)/(f1 - eps*f3);
	f4 = f(x4);
	if ( f4 == 0.0 ) return x4;
	// Contract the bracket.
	if ( f3*f2 < 0.0 ) {
	    if ( f4*f2 < 0.0 ) {
		x1 = x4; f1 = f4;
	    } else {
		x1 = x3; f1 = f3;
		x2 = x4; f2 = f4;
	    }
	} else {
	    if ( f4*f1 < 0.0 ) {
		x2 = x4; f2 = f4;
	    } else {
		x1 = x4; f1 = f4;
		x2 = x3; f2 = f3;
	    }
	}
    } // end while
    return x4;
} // end solve()


int bracket(std::function<double (double)>f, double &x1, double &x2,
	    int max_try, double factor)
{
    // Bracket a root of f(x).
    //
    // Input:
    //    f: user-supplied function f(x)
    //    x1: first end of range
    //    x2: other end of range
    // Returns:
    //    0: successfully bracketed a root
    //   -1: failed to bracket a root
    // On return (x1, x2) should bracketing the root.
    //
    if ( x1 == x2 ) {
	throw std::runtime_error("Bad initial range given to bracket.");
    }
    double f1 = f(x1);
    double f2 = f(x2);
    for ( int i = 0; i < max_try; ++i ) {
        if ( f1*f2 < 0.0 ) return 0; // we have success
        if ( fabs(f1) < fabs(f2) ) {
            x1 += factor * (x1 - x2);
            f1 = f(x1);
        } else {
            x2 += factor * (x2 - x1);
            f2 = f(x2);
	}
    }
    // If we leave the loop here, we were unsuccessful.
    return -1;
} // end bracket()
