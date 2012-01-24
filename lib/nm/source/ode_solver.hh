/*  \file ode_solver.hh
 *  \ingroup nm
 *  \brief Declarations for the generic ODE solver class.
 *  \author Rowan J Gollan
 *  \version 20-Feb-2006
 **/

#ifndef ODESOLVER_HH
#define ODESOLVER_HH

#include <iostream>
#include <string>
#include <valarray>
#include "no_fuss_linear_algebra.hh"
#include "ode_system.hh"
#include "ode_step.hh"


#define ODE_WARNINGS 1
#define MIN_H_FRACTION 1.0e-6

/** \brief Base class providing generic ode solving capabilities.

\author Rowan J Gollan
\version 20-Feb-2006

This base class provides the simplest set of functions
to solve an ODE over the interval \f$h_0\f$ --> \f$h_f\f$
based on a given stepping algorithm.
It accepts a number of arguments which allow for the tuning
of the solution algorithm.

**/

class OdeSolver {
public:
    /// \brief Default constructor
    OdeSolver();

    /// \brief Normal constructor
    OdeSolver( const std::string name, int ndim, const std::string step_name,
	       int max_step_attempts=4, double max_increase_factor=1.15,
	       double max_decrease_factor=1.0e-2, double decrease_factor=0.333 );

    /// \brief Copy constructor
    OdeSolver( const OdeSolver &o );

    /// \brief Default destructor
    virtual ~OdeSolver();

    /// \brief clone() function
    virtual OdeSolver* clone();

    void set_constants(const std::string name, int ndim, const std::string step_name,
		       int max_step_attempts, double max_increase_factor,
		       double max_decrease_factor, double decrease_factor,
		       double err_tol=1.0e-9);

    // -------- Services --------- //
    
    /// \brief string representation for the solver.
    virtual std::string str() const;

    /// \brief return the name of the solver.
    ///
    /// \author Rowan J Gollan
    /// \version 20-Feb-2006
    ///
    std::string name() const { return name_; }

    /// \brief return the dimensionality of the problem
    ///
    /// \author Rowan J Gollan
    /// \version 20-Feb-2006
    ///
    int ndim() const { return ndim_; }

    /// \brief return ther pointer to the ode step
    ///
    /// \author Rowan J Gollan
    /// \version 20-Feb-2006
    ///
    OdeStep* step() const { return step_; }

    /// \brief set the ode step type
    /// 
    /// \author Rowan J Gollan
    /// \version 21-Feb-2006
    /// 
    void set_step( OdeStep *step ) { delete step_; step_ = step->clone(); }

    // ------ Behaviour of an ODE solver ----- //

    /// \brief A procedure for solving the ODE over a specified interval.
    bool solve_over_interval( OdeSystem &ode, double x0, double xf, double *h,
			      const std::valarray<double> &yin, std::valarray<double> &yout );

protected:
    std::string name_;            ///< name for the solver
    int ndim_;                    ///< dimensionality of the problem
    int max_step_attempts_;       ///< number of attempts at taking a step.
    double max_increase_factor_;  ///< maximum allowable stepsize increase
    double max_decrease_factor_;  ///< maximum allowable stepzise decrease
    double decrease_factor_;      ///< decrease factor for failed step

    OdeStep *step_;               ///< a pointer to the stepping algorithm.
    
    std::valarray<double> y_save_; 
    std::valarray<double> y_work_;
};




#endif
