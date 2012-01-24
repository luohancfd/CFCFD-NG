/** \file zero_system.hh
 * \brief A class for creating a system for a zero finding technique.
 *
 * This object can be passed to a zero solver.  It encapsulates
 * the behaviour of a system of nonlinear functions, ie.
 * the function can be evaulated and so too can its Jacobian.
 *
 * \author Rowan J Gollan
 * \date 24-Apr-2006
 *
 **/

#ifndef ZERO_SYSTEM_HH
#define ZERO_SYSTEM_HH

#include <vector>
#include <valarray>

#include "no_fuss_linear_algebra.hh"

class ZeroSystem {
public:
    /// \brief Normal constructor
    ZeroSystem();

    /// \brief Copy constructor
    ZeroSystem( const ZeroSystem &z );

    /// \brief Default destructor
    virtual ~ZeroSystem();

    virtual int f( const std::valarray<double> &y, std::valarray<double> &G ) = 0;
    
    virtual int Jac( const std::valarray<double> &y, Valmatrix &dGdy ) = 0;

};

class ZeroFunction {
public:
    ZeroFunction();
    ZeroFunction(const ZeroFunction& zfun);
    virtual ~ZeroFunction();
    
    virtual int eval(double x, double &y) = 0; // y = f(x)
};

#endif
