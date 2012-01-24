/** \file quadrature.hh
 * \brief Declarations for various integration methods
 *
 * \author Daniel F Potter
 * \date 05-May-2009
 * \version 05-May-2009 - Just Gaussian 10 point quadrature for starters
 **/
 
#ifndef QUADRATURE_HH
#define QUADRATURE_HH
 
double gaussian_n10_integration( double (*f)( double x ), double a=-1.0, double b=1.0 );

#endif
