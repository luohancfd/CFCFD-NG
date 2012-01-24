/** \file quadrature.hh
 * \brief Definitions for various integration methods
 *
 * \author Daniel F Potter
 * \date 05-May-2009
 * \version 05-May-2009 - Just Gaussian 10 point quadrature for starters
 **/

using namespace std;
 
/// \brief Inspired by Numerical-Recipee's at http://www.fizyka.umk.pl/nrbook/c4-5.pdf p.148

double gaussian_n10_integration( double (*f)( double x ), double a, double b )
{
    // 0. Firstly define abscissas and weights (first entries are dummy values to make length even)
    static double x [] = { 0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285 };
    static double w [] = { 0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443 };
    
    // 1. Account for non [-1:1] range
    double xm = 0.5*(b+a);
    double xr = 0.5*(b-a);
    
    // 2. Perform integration, each mirrored pair at the same time
    double integral = 0.0;
    double dx;
    for ( int j=1; j<=5; ++j ) {
	dx = xr*x[j];
	integral += w[j]*( (*f)(xm+dx) + (*f)(xm-dx) );
    }
    
    // 3. Scale the answer to the range of integration
    return integral *= xr;
}
