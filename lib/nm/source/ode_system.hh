/*  \file ode_system.hh
 *  \ingroup nm
 *  \brief Declarations for an ode system used in an ode solver.
 *  \author Rowan J Gollan
 *  \version 20-Feb-2006
 */

#ifndef ODESYS_HH
#define ODESYS_HH

#include <iostream>
#include <string>
#include <vector>
#include "no_fuss_linear_algebra.hh"

/** \brief Base (abstract) class of an ode system.

\author Rowan J Gollan
\version 20-Feb-2006 : initial coding

This abstract class encapsulates the behaviour expected
of an ODE system.
Along with providing a function evaluation for the ODE
system, it is desirable that the system
has some means of estimating an appropriate timestep.

Also, in some instances, the system can check that the answer
is sensible at the end of the ODE stepping routine.
For example, in a finite-rate chemistry ODE problem, mass
should be conserved at the end of ODE stepping.

The ODE system can also hold some information about various options
that can be used to solve the system.

**/

class OdeSystem {
public:
    /// \brief Default constructor
    OdeSystem();

    /// \brief Normal constructor
    OdeSystem( int ndim, bool system_test  );

    /// \brief Copy constructor
    OdeSystem( const OdeSystem &o );

    /// \brief Default destructor
    virtual ~OdeSystem();

    void set_constants(int ndim, bool system_test);

    // ---- Services ---- //
    
    /// \brief string representation of the system
    virtual std::string str() const;

    /// \brief return boolean _apply_system_test
    /// 
    /// \author Rowan J Gollan
    /// \version 20-Feb-2006
    ///
    bool apply_system_test() { return apply_system_test_; }

    /// \brief return system dimensionality
    ///
    /// \author Rowan J Gollan
    /// \version 20-Feb-2006
    ///
    int ndim() { return ndim_; }

    // --- Functions defining the behaviour --- //

    /// \brief The RHS evaluation for the ODE
    /// 
    /// \author Rowan J Gollan
    /// \version 20-Feb-2006
    ///
    /// This function represents the basic functionality
    /// of an ODE system.  Given the generic ODE system,
    /// \f[ \frac{dev\vec{y}}{dt} = f(\vec{y}) \f],
    /// then this eval returns the result of \f$f(\vec{y})\f$.
    ///
    virtual int eval( const std::vector<double> &y, std::vector<double> &ydot ) = 0;

    /// \brief A RHS evaluation which forward and backward production
    virtual int eval_split( const std::vector<double> &y, 
			    std::vector<double> &q, std::vector<double> &L );

    /// \brief A stepsize selection function.
    virtual double stepsize_select( const std::vector<double> &y );

    /// \brief A check for the accuracy of the result
    ///        based on system constraints.
    virtual bool passes_system_test( std::vector<double> &y );

    /// \brief A flag indicating if the RHS has been evaluated.
    ///
    /// Some complex ODE systems don't need to re-evaluate 
    /// all of the internal data at every timestep.
    /// This public flag is provided so that the stepping
    /// mechanism can set and retrieve this value.
    bool called_at_least_once;

protected:
    int ndim_;               ///< dimensionality of ODE system
    bool apply_system_test_; ///< indicates if a system test should be applied
};

#ifndef SWIG
/// \brief Overloaded ostream operator
std::ostream& operator<<( std::ostream &os, const OdeSystem &o );
#endif

#endif // ODESYS_HH
