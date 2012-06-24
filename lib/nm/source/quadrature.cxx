/** \file quadrature.hh
 * \brief Definitions for various integration methods
 *
 * \author Daniel F Potter
 * \date 05-May-2009
 * \version 05-May-2009 - Just Gaussian 10 point quadrature for starters
 **/

#include <math.h>
#include <iostream>

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

/// \brief Shamefully borrowed from wikipedia: http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method

double adaptiveSimpsonsAux(double (*f)(double x, void * params), double a, double b, void * params, double epsilon,
                         double S, double fa, double fb, double fc, int bottom)
{
  cout << "bottom = " << bottom << endl;
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
  double fd = f(d,params), fe = f(e,params);
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
    return S2 + (S2 - S)/15;
  return adaptiveSimpsonsAux(f, a, c, params, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
         adaptiveSimpsonsAux(f, c, b, params, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

//
// Adaptive Simpson's Rule
//
double adaptiveSimpsons(double (*f)(double x, void * params),   // ptr to function
                           double a, double b,  // interval [a,b]
                           void * params,   // parameters for the derivative function
                           double epsilon,  // error tolerance
                           int maxRecursionDepth) {   // recursion cap
  double c = (a + b)/2, h = b - a;
  double fa = f(a,params), fb = f(b,params), fc = f(c,params);
  double S = (h/6)*(fa + 4*fc + fb);
  return adaptiveSimpsonsAux(f, a, b, params, epsilon, S, fa, fb, fc, maxRecursionDepth);
}
