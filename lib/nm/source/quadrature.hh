/** \file quadrature.hh
 * \brief Declarations for various integration methods
 *
 * \author Daniel F Potter
 * \date 05-May-2009
 * \version 05-May-2009 - Just Gaussian 10 point quadrature for starters
 *          11-Apr-2012 - Adaptive Simpson's rule
 **/
 
#ifndef QUADRATURE_HH
#define QUADRATURE_HH
 
double gaussian_n10_integration( double (*f)( double x ), double a=-1.0, double b=1.0 );

double adaptiveSimpsonsAux(double (*f)(double x, void * params), double a, double b,
                           void * params, double epsilon, double S, double fa,
                           double fb, double fc, int bottom);

double adaptiveSimpsons(double (*f)(double x, void * params),   // ptr to function
                           double a, double b,  // interval [a,b]
                           void * params,   // parameters for the derivative function
                           double epsilon,  // error tolerance
                           int maxRecursionDepth);

#endif
