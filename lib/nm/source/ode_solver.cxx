/*  \file ode_solver.cxx
 *  \ingroup nm
 *  \brief Definitions for the generic ODE solver class.
 *  \author Rowan J Gollan
 *  \version 20-Feb-2006
 **/
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include "no_fuss_linear_algebra.hh"
#include "ode_system.hh"
#include "ode_step.hh"
#include "ode_solver.hh"
#include "../../util/source/useful.h"

using namespace std;

/// \brief Default constructor
///
/// \author Rowan J. Gollan
/// \version 07-Jan-07
///
OdeSolver::OdeSolver() {}

/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
/// \param name      : a name for the ode solver
/// \param ndim      : the dimensionality for the problem
/// \param step_name : the desired stepping algorithm
///
OdeSolver::OdeSolver( const string name, int ndim, const string step_name,
		      int max_step_attempts, double max_increase_factor,
		      double max_decrease_factor, double decrease_factor )
    : name_( name ), ndim_( ndim ), max_step_attempts_( max_step_attempts ),
      max_increase_factor_( max_increase_factor ), max_decrease_factor_( max_decrease_factor ),
      decrease_factor_( decrease_factor)
{

    if( step_name == "euler" ) {
	step_ = new EulerStep("euler", ndim);
    }
    else if( step_name == "qss" ) {
	step_ = new QssStep("qss", ndim, 10, 1.0e-5, 3.0 );
    }
    else if( step_name == "rkf" ) {
	step_ = new RKFStep( "rkf", ndim, 1.0e-9);
    }
    else if ( step_name == "dp853" ) {
	step_ = new DP853Step("dp853", ndim, 1.0e-9 );
    }
    else {
	cout << "OdeSolver::OdeSolver() : problem during initialisation.\n"
	     << "The step_name= " << step_name
	     << " is unknown or not yet implemented.\n"
	     << "Bailing out.\n";
	exit(BAD_INPUT_ERROR);
    }

    y_save_.resize(ndim);
    y_work_.resize(ndim);

}

/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
OdeSolver::OdeSolver( const OdeSolver &o )
    : name_( o.name_ ), ndim_( o.ndim_ ),
      max_step_attempts_( o.max_step_attempts_ ), max_increase_factor_(o.max_increase_factor_),
      max_decrease_factor_( o.max_decrease_factor_ ), decrease_factor_( o.decrease_factor_)
     
{
    step_ = o.step()->clone();
    y_save_.resize(ndim_);
    y_work_.resize(ndim_);
}

/// \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
/// We need to delete the memory allocated for step_.
///
OdeSolver::~OdeSolver()
{
    delete step_;
}

/// \brief clone() function
///
/// \author Rowan J Gollan
/// \version 25-Feb-2006
///
OdeSolver* OdeSolver::clone()
{
    return new OdeSolver(*this);
}

/// \brief Set constants for ode solver.
/// 
/// \author Rowan J. Gollan
/// \version 07-Jan-07
///
void OdeSolver::set_constants(const string name, int ndim, const string step_name,
			      int max_step_attempts, double max_increase_factor,
			      double max_decrease_factor, double decrease_factor,
			      double err_tol)
{
    name_ = name;
    ndim_ = ndim;
    max_step_attempts_ =  max_step_attempts;
    max_increase_factor_ = max_increase_factor;
    max_decrease_factor_ = max_decrease_factor;
    decrease_factor_ = decrease_factor;

    if( step_name == "euler" ) {
	step_ = new EulerStep("euler", ndim);
    }
    else if( step_name == "qss" ) {
	step_ = new QssStep("qss", ndim, 10, err_tol, 3.0 );
    }
    else if( step_name == "rkf" ) {
	step_ = new RKFStep( "rkf", ndim, err_tol);
    }
    else if ( step_name == "dp853" ) {
	step_ = new DP853Step("dp853", ndim, err_tol);
    }
    else {
	cout << "OdeSolver::set_constants() : problem during setting of constants.\n"
	     << "The step_name= " << step_name
	     << " is unknown or not yet implemented.\n"
	     << "Bailing out.\n";
	exit(BAD_INPUT_ERROR);
    }

    y_save_.resize(ndim);
    y_work_.resize(ndim);
}

/// \brief string representation for the solver.
/// 
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
string OdeSolver::str() const
{
    ostringstream ost;
    ost << "OdeSolver(\n"
	<< "   name= " << name_ << endl
	<< "   ndim= " << ndim_ << endl
        << "   step= " << step_->name() << endl
	<< "   max_step_attempts= " << max_step_attempts_ << endl
	<< "   max_increase_factor= " << max_increase_factor_ << endl
	<< "   max_decrease_factor= " << max_decrease_factor_ << endl
	<< "   decrease_factor= " << decrease_factor_ << endl
	<< ")\n";
    return ost.str();
}

/// \brief A procedure for solving the ODE over a specified interval.
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
/// \param ode   : the ode system being solved (provides RHS eval)
/// \param x0    : the initial starting point for the independent variable
///                (I've tried to refrain from calling this t0 as there
///                 are ODEs that can be solved spatially also.)
/// \param xf    : the final point for solution termination.
/// \param h     : if <= 0.0, a step size selection will be made, otherwise
///                the value in h is used as the initial step size.
///                At the end of the function, this holds the last used
///                timestep (which may be useful for later restart)
/// \param yin   : the input vector
/// \param yout  : the output vector
///
/// \return      : true if successful, otherwise false.
///                If false, the value in yout should not be used.
///
/// The general ODE solving procedure will take multiple internal
/// steps to reach the desired solution value.
///
bool OdeSolver::solve_over_interval( OdeSystem &ode, double x0, double xf, double *h,
				     const valarray<double> &yin, valarray<double> &yout )
{
    // 1. Check the bounds.
    double hdiff = xf - x0;

    if( hdiff <= 0.0 ) {
	cout << "OdeSolver::solve_over_interval() \n"
	     << "Warning - the interval " << x0 << " to " << xf 
	     << " is not positive, hdiff = " << hdiff << endl
	     << "Bailing out.\n";
	exit(BAD_INPUT_ERROR);
    }

    // 2. Select a step size
    if( *h <= 0.0 ) { // We need to select a step size
	*h = ode.stepsize_select( yin );
    }

    if( *h > hdiff ) {
	*h = hdiff; // No point stepping large than solution interval.
    }

    // 3. Keep a copy of the input in case we need to restart
    copy_valarray( yin, y_save_ );
    copy_valarray( yin, y_work_ );


    // 4. Attempt a solution
    double xcurr = x0;
    double h_last;
    double h_final = *h;
    int steps = 0;

    const int max_system_failures = 5;
    int system_failures = 0;

    while( *h < 0.5*(xf - xcurr)) {
	int step_attempt = 0;
	//cout << "h= " << *h << endl;
	for( ; step_attempt < max_step_attempts_; ++step_attempt ) {
	    
	    //cout << "step_attempt= " << step_attempt << endl;
	    h_last = *h;
	    
	    if( step_->advance( ode, y_work_, yout, h ) ) { // step is successful
		steps++;
		//cout << "Step successful.\n";
		if( ode.apply_system_test() && !ode.passes_system_test( yout ) ) {
		    //cout << "BUT DID NOT PASS SYSTEM TEST.\n";
		    system_failures++;
		    //cout << "system_failures= " << system_failures << endl;
		    if ( system_failures == max_system_failures )
			return false;
		    // tolerance on step not sufficient
		    *h = h_last * decrease_factor_;
		    // Often the step will fail irrespective of the step size,
		    // so we only want to attempt a finite number of reductions
		    // --step_attempt;
		    // --steps;
		    //cout << "Now h= " << *h << endl;
		    continue;
		}
		// this is what was replaced... just in case.
		// if( ode.apply_system_test() ) {
		//     if( ode.passes_system_test( yout ) ) {
		// 	;
		//     }
		//     else {
		// 	*h = h_last * decrease_factor_;
		// 	continue;
		//     }
		// }
		xcurr += h_last;

		copy_valarray( yout, y_work_ );
		
		h_final = h_last;
		
		// Now make sure the next step is reasonable...
		if( *h < h_last ) {
		    // We just had a successful step so it should be safe to
		    // stay the same
		    *h = h_last;
		}
		else if( *h > ( max_increase_factor_ * h_last ) ) {
		    // We don't want to increase too wildly.
		    *h = max_increase_factor_ * h_last;
		}
		break;
	    }
	    else { // The step failed.
		//cout << "STEP FAILED.\n";
		// So we select a new timestep and start again...
		if( *h >= h_last ) {
		    // For some reason the suggested timestep is the same or larger.
		    *h = h_last * decrease_factor_;
		}
		else if( *h < (h_last * max_decrease_factor_) ) {
		    // The reduced step is going to be too small.
		    *h = h_last * max_decrease_factor_;
		}
		else if ( !finite(*h) ) {
		    // The time-step is not sensible
		    *h = h_last * max_decrease_factor_;
		}
		// impose MIN_H_FRACTION constraint on *h
		if ( *h / hdiff < MIN_H_FRACTION ) {
		    // The reduced step is going to be too small
		    *h = hdiff * MIN_H_FRACTION;
		}
		//cout << "So next step h= " << *h << endl;
	    }
	}

	if( step_attempt >= max_step_attempts_ ) {
#           if ODE_WARNINGS >= 2	    
	    cout << "The ODE solver failed trying to take a single step\n"
		 << "after " << max_step_attempts_ << " attempts." << endl;
#           endif
	    return false;
	}
	
    }

    // Now we have two steps left to take
    // So we take two half steps.
    h_final = *h;
    double h_trial = (xf - xcurr);
    int no_steps;
    int step_attempt = 0;
    for( ; step_attempt < max_step_attempts_; ++step_attempt ) {
	no_steps = int(2*pow(2.0, step_attempt));
	*h = h_trial / no_steps;
	int i = 0;
	for( ; i < no_steps; ++i ) {
	    if( step_->advance( ode, y_work_, yout, h ) ) {
		if( steps == 0 ) { // We've never had a h suggestion.
		    h_final = *h;
		}
		copy_valarray( yout, y_work_ );
		*h = h_trial / no_steps;
	    }
	    else
		break;
	}
	if( i == no_steps ) {
	    // Steps ok, do system check...
	    if( ode.apply_system_test() && !ode.passes_system_test( yout ) ) {
		//cout << "OdeSolver::solve_over_interval():" << endl
		//     << "    system check failed at the end of the interval." << endl
		//     << "    Maybe we should be exiting to the operating system." << endl;
		return false;
		// exit(1);
	    }
	    else
		break;
	}
    }

    if( step_attempt >= max_step_attempts_ ) {
#       if ODE_WARNINGS >= 2
	cout << "The ODE solver failed on the finishing steps.\n"
	     << "after " << max_step_attempts_ << " attempts." << endl;
#       endif
	return false;
    }
	
    // 5. There should be no need for interpolation as we always take
    //    a step to bring us to the finishing value.
  
    // If we've made it this far we must be successful.
    *h = h_final;

    return true;
}

