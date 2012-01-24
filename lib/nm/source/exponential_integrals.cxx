// Author: Daniel F. Potter
// Date: 28-Apr-2009

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../util/source/useful.h"
#include "exponential_integrals.hh"

using namespace std;

/* Exponential-integration approximations */

// NOTE: E_n() is useful for testing purposes, but must test for 'n'

double E_n( int n, double x )
{
    /* Curve-fits to nth order exponential integrals of the form:
       E_n( x ) = int [1, inf] ( w^(-n) exp( -w x ) dw
       Obtained from the following paper by Chris Johnston:
       "Radiative Heating Methodology for the Huygens Probe" AIAA 2006-3426 */
       
    if      (n==1) return 2.5910 * exp( - 18.700 * x ) + 1.7080 * exp( -2.110 * x );
    else if (n==2) return 0.2653 * exp( -  8.659 * x ) + 0.7347 * exp( -1.624 * x );
    else if (n==3) return 0.0929 * exp( -  4.080 * x ) + 0.4071 * exp( -1.330 * x );
    else {
	cout << "The requested exponential integral is not implemented" << endl;
	exit( FAILURE );
    }
}

// NOTE: explicit functions for numerical applications

double E_1( double x )
{
    // Johnston (2006) approximation
    return 2.5910 * exp( - 18.700 * x ) + 1.7080 * exp( -2.110 * x );

}

double E_2( double x )
{
    // Johnston (2006) approximation
    return 0.2653 * exp( -  8.659 * x ) + 0.7347 * exp( -1.624 * x );

}

# define E3_METHOD 0

double E_3( double x )
{
#   if E3_METHOD == 0
    // Johnston (2006) approximation
    return 0.0929 * exp( -  4.080 * x ) + 0.4071 * exp( -1.330 * x );
#   elif E3_METHOD == 1
    // Karl (2001) approximation 
    if ( x == 0.0 ) return 0.5;
    else return 0.5 - x + 0.5 * ( 0.9228 - log(x) ) * x * x + x * x * x / ( 6.0 );
#   endif
}

# undef E3_METHOD

