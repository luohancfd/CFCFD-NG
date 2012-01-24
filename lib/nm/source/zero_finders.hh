/** \file zero_finders.hh
 * \ingroup nm
 * \brief A collection of zero-finding methods for multi-dimensional
 *        systems.
 *
 * This module provides a collection of zero-finding methods for solving
 * systems of nonlinear equations.
 *
 * \date 21-Apr-2006
 * \author Rowan J Gollan
 * \version 21-Apr-2006 - initial coding
 *
 **/

#ifndef ZERO_FINDERS_HH
#define ZERO_FINDERS_HH

#include <vector>
#include <valarray>

#include "no_fuss_linear_algebra.hh"
#include "zero_system.hh"

/** \brief Base clase defining generic zero_finder behaviours.

\author Rowan J Gollan
\version 21-Apr-2006

This base class encapsulates the minimal set of behviours
expected of a zero-finder.

**/

class ZeroFinder {
public:
    /// \brief Default constructor
    ZeroFinder();

    /// \brief Normal constructor
    ZeroFinder( int ndim, double tolerance, int max_iterations );

    /// \brief Copy constructor
    ZeroFinder( const ZeroFinder &z );

    /// \brief Default destuctor
    virtual ~ZeroFinder();

    /// \brief Set constants for the zero finder.
    void set_constants(int ndim, double tolerance, int max_iterations);

    //---- Defining behaviour ----//

    /// \brief Solve the system
    virtual int solve( ZeroSystem &zsys, const std::valarray<double> &y_guess, std::valarray<double> &y_out ) = 0;

protected:
    int ndim_;     ///< Number of variables to solve for
    double tol_;   ///< acceptable tolerance for stopping
    int max_iter_; ///< maximum number of iterations while attempting solve
};

/// \brief  Bisection method
/// \author Brendan T. O'Flaherty 

class Bisection {
public:
    Bisection();
    Bisection(ZeroFunction *zfun, double tol);
    Bisection(const Bisection &b);
    virtual ~Bisection();

    Bisection* clone() const;

    // Bisection method is executed by providing bounds.
    double operator()(double x_left, double x_right);
    
    // change the function later
    void set_zfun(ZeroFunction *zfun) { zfun_ = zfun; }

private:
    double tol_;
    ZeroFunction *zfun_;
};

// \brief Muller's method
// \author Brendan T. O'Flaherty

class Muller {
public:
    Muller();
    Muller(ZeroFunction *zfun, double tol);
    Muller(const Muller &n);
    virtual ~Muller();
    
    Muller* clone() const;
    
    // Muller method is executed by providing bounds.
    double operator()(double x_left, double x_right);
    
    // change the function later
    void set_zfun(ZeroFunction *zfun) { zfun_ = zfun; }
    
private:
    double tol_;
    double x_mid_, y_mid_;
    ZeroFunction *zfun_;
};

class NewtonRaphsonZF : public ZeroFinder {
public:
    /// \brief Default constructor
    NewtonRaphsonZF();
    
    /// \brief Normal constructor
    NewtonRaphsonZF( int ndim, double tolerance, int max_iterations, bool has_Jacobian=true );

    /// \brief Copy constructor
    NewtonRaphsonZF( const NewtonRaphsonZF &n );

    /// \brief Default destructor
    virtual ~NewtonRaphsonZF();

    /// \brief Set constants for zero finder.
    void set_constants(int ndim, double tolerance, int max_iterations, bool has_Jacobian=true);

    //---- Defining behaviour ----//
    int solve( ZeroSystem &zsys, const std::valarray<double> &y_guess, std::valarray<double> &y_out );

private:
    bool has_Jac_;
    std::valarray<double> G_;
    std::valarray<double> minusG_;
    Valmatrix dGdy_;
    std::valarray<double> y_old_;
    std::valarray<double> y_new_;
    std::valarray<double> dely_;

    bool test_tol();
};

#endif
  
