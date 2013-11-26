/** \file zero_system.hh
 * \brief A class for creating a system for a zero finding technique.
 *
 * This object can be passed to a zero solver.  It encapsulates
 * the behaviour of a system of nonlinear functions, ie.
 * the function can be evaulated and so too can its Jacobian.
 *
 * \author Rowan J Gollan
 * \date 24-Apr-2006
 * \version 26-Nov-2013 -- change to using vector (PJ)
 *
 **/

#ifndef ZERO_SYSTEM_HH
#define ZERO_SYSTEM_HH

#include <vector>

#include "no_fuss_linear_algebra.hh"

class ZeroSystem {
public:
    /// \brief Normal constructor
    ZeroSystem();

    /// \brief Copy constructor
    ZeroSystem( const ZeroSystem &z );

    /// \brief Default destructor
    virtual ~ZeroSystem();

    virtual int f( const std::vector<double> &y, std::vector<double> &G ) = 0;
    
    virtual int Jac( const std::vector<double> &y, Valmatrix &dGdy ) = 0;

};

class ZeroFunction {
public:
    ZeroFunction();
    ZeroFunction(const ZeroFunction& zfun);
    virtual ~ZeroFunction();
    
    virtual int eval(double x, double &y) = 0; // y = f(x)
};

#endif
