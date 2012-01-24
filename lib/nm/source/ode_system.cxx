/*  \file ode_system.cxx
 *  \ingroup nm
 *  \brief Definitions for an ode system used in an ode solver.
 *  \author Rowan J Gollan
 *  \version 20-Feb-2006
 */

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "no_fuss_linear_algebra.hh"
#include "ode_system.hh"
using namespace std;

/// \brief Default constructor
/// 
/// \author Rowan J. Gollan
/// \version 07-Jan-07
///
OdeSystem::OdeSystem() {}


/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
/// \param ndim        : dimensionality of ODE system
/// \param system_test : flag indicating if a system test is used
///
OdeSystem::OdeSystem( int ndim, bool system_test )
    : ndim_( ndim ), apply_system_test_( system_test ) {}

/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
OdeSystem::OdeSystem( const OdeSystem &o )
    : ndim_( o.ndim_ ), apply_system_test_( o.apply_system_test_ ) {}

/// \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
OdeSystem::~OdeSystem() {}

/// \brief Set constants for the ode system.
///
/// \author Rowan J. Gollan
/// \version 07-Jan-07
///
void OdeSystem::set_constants(int ndim, bool system_test)
{
    ndim_ = ndim;
    apply_system_test_ = system_test;
}

/// \brief string representation of the system
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
string OdeSystem::str() const
{
    ostringstream ost;
    ost << "OdeSystem(\n"
	<< "   ndim= " << ndim_ << endl
	<< "   apply_system_test= " << apply_system_test_ << endl
	<< ")\n";

    return ost.str();
}

/// \brief A RHS evaluation which forward and backward production
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
/// 
/// Some ODE integrators make use of the positive (production)
/// and negative (loss) terms in their stepping algorithms.
/// This evaulation function computes those behaviours.
///
/// If the ODE system can't provide those values, we can
/// use this generic function to find y' and then
/// assign it as a production or loss depending on its sign.
///
int OdeSystem::eval_split( const valarray<double> &y,
			   valarray<double> &q, valarray<double> &L )
{
    // Use q to hold values rather than grabbing memory.
    eval(y, q);
    for( int i = 0; i < ndim_; ++i ) {
	if( q[i] >= 0.0 ) {
	    L[i] = 0.0;
	}
	else {
	    L[i] = -q[i];
	    q[i] = 0.0;
	}
    }

    return 0;

}

/// \brief A stepsize selection function.
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
/// Trying to define a function that will select a stable
/// step size for an ODE system is crazy.
/// The step size selection depends both on the system
/// and the integration method that is being used.
/// This trivial function just sets the step size to
/// 1.0e-6.  This will be completely inappropriate
/// for 99.99% of cases that are passed to the
/// ODE solver.
///
/// \note The implementer of derived classes of OdeSystem
///       is strongly encouraged to code their own
///       step size selection function depending on their
///       ODE system.
///
double OdeSystem::stepsize_select( const valarray<double> &y )
{
    const double silly_number_for_timestep = 1.0e-6;
    return silly_number_for_timestep;
}

/// \brief A check for the accuracy of the result
///        based on system constraints.
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
/// This implementation does nothing but just returns
/// true.  It is provided so that the derived class implementer
/// is not forced to implement a passes_system_test() function
/// if it's not appropriate for his/her system.
///
bool OdeSystem::passes_system_test( valarray<double> &y )
{
    return true;
}

/// \brief Overloaded ostream operator
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
ostream& operator<<( ostream &os, const OdeSystem &o )
{
    os << o.str();
    return os;
}
