/*  \file ode_step.hh
 *  \ingroup nm
 *  \brief Declarations for the OdeStep class.
 *  \author Rowan J Gollan
 *  \version 20-Feb-2006
 **/

#ifndef ODESTEP_HH
#define ODESTEP_HH

#include <iostream>
#include <string>
#include <valarray>
#include "no_fuss_linear_algebra.hh"
#include "ode_system.hh"

/** \brief Base (abstract) class for an ode-stepping algorithm.

\author Rowan J Gollan
\version 20-Feb-2006 : initial coding.

This abstract coding encapsulates the behaviour of an
ODE stepping algorithm.
The minimal functionality provided increments and ODE system
a certain amount, dt.

Additionally, some stepping algorithms can report on the 
accuracy of their stepping (ie. if it was successful or
not) and suggest a new timestep for the next step.
This functionality is also encapsulated.

**/

class OdeStep {
public:
    /// \brief Normal constructor
    OdeStep( const std::string name, int ndim );

    /// \brief Copy constructor
    OdeStep( const OdeStep &o );
    
    /// \brief Default destructor
    virtual ~OdeStep();

    // ------- Services ---------- //
    
    /// \brief string representation of the step
    virtual std::string str() const;

    /// \brief allocate new memory and make a clone
    ///        of an ode step
    ///
    /// \author Rowan J Gollan
    /// \version 20-Feb-2006
    ///
    /// This is a pure virtual function.  Derived classes must provide
    /// this function.
    ///
    virtual OdeStep* clone() = 0;

    /// \brief return the name of the alogorithm
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
    int ndim() const { return ndim_; };

    // ------ Behaviour of an ODE step ----- //
    
    /// \brief Takes an incremental ODE step
    /// 
    /// \author Rowan J Gollan
    /// \version 20-Feb-2006
    ///
    /// \param ode  : an OdeSystem which encapsualtes auxilliary data
    ///               and behaviour of the system.
    /// \param yin  : the vector of doubles containing current input
    /// \param yout : this vector has the values of the output at the
    ///               the end of the step.
    /// \param h    : the size of the step - this is a pointer as it will
    ///               be updated with a suggested step size for the next step.
    ///
    /// \returns    : true if step is numerically successful, false otherwise
    ///               if false, we shouldn't use the values in yout
    ///
    /// \note This is a pure virtual function.  The abstract class defines
    ///       how this function should behave but gives no particular implementation.
    ///
    virtual bool advance( OdeSystem &ode, const std::valarray<double> &yin,
			  std::valarray<double> &yout, double *h ) = 0;

protected:
    std::string name_;   ///< name for the stepping algorithm
    int ndim_;           ///< dimensionality of the problem
};

/// \brief Overloaded ostream operator
std::ostream& operator<<( std::ostream &os, const OdeStep &o );

//-----------------------------------------------------------------------//
// Derived OdeSteps: actual ODE stepping algorithms                      //
//-----------------------------------------------------------------------//

/** \brief Declarations for the EulerStep class.

\author Rowan J Gollan
\version 20-Feb-2006

This class implements the time-honoured Euler update for
an ODE system.
\f[ y_{n+1} = y_n + h y' \f]

**/

class EulerStep : public OdeStep {
public:
    /// \brief Normal constructor
    EulerStep( const std::string name, int ndim );

    /// \brief Copy constructor
    EulerStep( const EulerStep &o );
    
    /// \brief Default destructor
    virtual ~EulerStep();

    /// \brief Clone of Euler step.
    EulerStep* clone();

    // ------ Behaviour of an Euler step ----- //
    
    /// \brief Takes an Euler step.
    bool advance( OdeSystem &ode, const std::valarray<double> &yin,
		  std::valarray<double> &yout, double *h );
private:
    std::valarray<double> y_dot_;
};

/** \brief Declarations for the ModEulerStep class.

\author Rowan J Gollan
\version 21-Feb-2006

This class implements the modified Euler method (sometimes
known as the improved Euler method).
The implementation here is an explicit formulation
and the update is computed as:

\f[ y_{n+1} = y_n + h \frac{y'_n + y'_{n+1}}{2} \f]

where \f$y'_{n+1}\f$ is evaluated in an explicit manner, ie:

\f[ y'_{n+1} = f(y_{n+1}) = f( y_n + h f(y_n) ) \f]

**/

class ModEulerStep : public OdeStep {
public:
    /// \brief Normal constructor
    ModEulerStep( const std::string name, int ndim );

    /// \brief Copy constructor
    ModEulerStep( const ModEulerStep &m );

    /// \brief Default destructor
    virtual ~ModEulerStep();

    /// \brief Clone of the ModEuler object
    ModEulerStep* clone();

    // -------- Behaviour of a modified Euler step ------ //

    /// \brief Takes a modified Euler step
    bool advance( OdeSystem &ode, const std::valarray<double> &yin, 
		  std::valarray<double> &yout, double *h );

private:
    std::valarray<double> dy_n_;
    std::valarray<double> dy_nplus1_;
    std::valarray<double> y_nplus1_;

};

/** \brief Declarations for the Runge-Kutta-Fehlberg algorithm

\author Rowan J Gollan
\version 21-Feb-2006

This class implements the Runge-Kutta-Fehlberg algorithm
for ODE integration.
This method can give an estimated stable step size
based on acceptable error.

**/

class RKFStep : public OdeStep {
public:
    /// \brief Normal constructor
    RKFStep( const std::string name, int ndim, double tol=1.0e-8 );
    
    /// \brief Copy constructor
    RKFStep( const RKFStep &r );

    /// \brief Default destructor
    virtual ~RKFStep();

    /// \brief Clone of the ModEuler object
    RKFStep* clone();

    // -------- Behaviour of a modified Euler step ------ //

    /// \brief Takes a step using the Runge-Kutta-Fehlberg method
    bool advance( OdeSystem &ode, const std::valarray<double> &yin, 
		  std::valarray<double> &yout, double *h );

private:
    void set_constants();
    double tol_;
    // Some working arrays for the model
    std::valarray<double> tmp_;
    std::valarray<double> yerr_;
    std::valarray<double> k1_;
    std::valarray<double> k2_;
    std::valarray<double> k3_;
    std::valarray<double> k4_;
    std::valarray<double> k5_;
    std::valarray<double> k6_;

    // Constants from
    // Cash and Karp (1990)
    // A Variable Order Runge-Kutta Method for
    // Initial Value Problems woth Rapidly Varying
    // Right-Hand Sides
    // ACM Transactions on Mathematical Software, 16:3 pp 201--222
    
    double a21_;
    double a31_;
    double a32_;
    double a41_;
    double a42_;
    double a43_;
    double a51_;
    double a52_;
    double a53_;
    double a54_;
    double a61_;
    double a62_;
    double a63_;
    double a64_;
    double a65_;

    // For autonomous systems (which we treat here)
    // these are not needed.
//     const double c2_ = 1.0/5.0;
//     const double c3_ = 3.0/10.0;
//     const double c4_ = 3.0/5.0;
//     const double c5_ = 1.0;
//     const double c6_ = 7.0/8.0;

    double b51_;
    double b53_;
    double b54_;
    double b56_;
    
    double b41_;
    double b43_;
    double b44_;
    double b45_;
    double b46_;

};

class DP853Step : public OdeStep {
public:
    /// \brief Normal constructor
    DP853Step( const std::string name, int ndim, double tol=1.0e-8 );
    
    /// \brief Copy constructor
    DP853Step( const DP853Step &d );

    /// \brief Default destructor
    virtual ~DP853Step();

    /// \brief Clone of the ModEuler object
    DP853Step* clone();

    // -------- Behaviour of a modified Euler step ------ //

    /// \brief Takes a step using the Runge-Kutta-Fehlberg method
    bool advance( OdeSystem &ode, const std::valarray<double> &yin, 
		  std::valarray<double> &yout, double *h );

private:
    void set_constants();
    double tol_;
    // Some working arrays for the model
    std::valarray<double> tmp_;
    std::valarray<double> yerr_;
    std::valarray<double> yerr2_;
    std::valarray<double> k1_;
    std::valarray<double> k2_;
    std::valarray<double> k3_;
    std::valarray<double> k4_;
    std::valarray<double> k5_;
    std::valarray<double> k6_;
    std::valarray<double> k7_;
    std::valarray<double> k8_;
    std::valarray<double> k9_;
    std::valarray<double> k10_;

    double b1_;
    double b6_;
    double b7_;
    double b8_;
    double b9_;
    double b10_;
    double b11_;
    double b12_;

    double bhh1_;
    double bhh2_;
    double bhh3_;

    double er1_;
    double er6_;
    double er7_;
    double er8_;
    double er9_;
    double er10_;
    double er11_;
    double er12_;

    double a21_;
    double a31_;
    double a32_;
    double a41_;
    double a43_;
    double a51_;
    double a53_;
    double a54_;
    double a61_;
    double a64_;
    double a65_;
    double a71_;
    double a74_;
    double a75_;
    double a76_;
    double a81_;
    double a84_;
    double a85_;
    double a86_;
    double a87_;
    double a91_;
    double a94_;
    double a95_;
    double a96_;
    double a97_;
    double a98_;
    double a101_;
    double a104_;
    double a105_;
    double a106_;
    double a107_;
    double a108_;
    double a109_;
    double a111_;
    double a114_;
    double a115_;
    double a116_;
    double a117_;
    double a118_;
    double a119_;
    double a1110_;
    double a121_;
    double a124_;
    double a125_;
    double a126_;
    double a127_;
    double a128_;
    double a129_;
    double a1210_;
    double a1211_;
    double a141_;
    double a147_;
    double a148_;
    double a149_;
    double a1410_;
    double a1411_;
    double a1412_;
    double a1413_;
    double a151_;
    double a156_;
    double a157_;
    double a158_;
    double a1511_;
    double a1512_;
    double a1513_;
    double a1514_;
    double a161_;
    double a166_;
    double a167_;
    double a168_;
    double a169_;
    double a1613_;
    double a1614_;
    double a1615_;

    double d41_;
    double d46_;
    double d47_;
    double d48_;
    double d49_;
    double d410_;
    double d411_;
    double d412_;
    double d413_;
    double d414_;
    double d415_;
    double d416_;
    double d51_;
    double d56_;
    double d57_;
    double d58_;
    double d59_;
    double d510_;
    double d511_;
    double d512_;
    double d513_;
    double d514_;
    double d515_;
    double d516_;
    double d61_;
    double d66_;
    double d67_;
    double d68_;
    double d69_;
    double d610_;
    double d611_;
    double d612_;
    double d613_;
    double d614_;
    double d615_;
    double d616_;
    double d71_;
    double d76_;
    double d77_;
    double d78_;
    double d79_;
    double d710_;
    double d711_;
    double d712_;
    double d713_;
    double d714_;
    double d715_;
    double d716_;
    
};

/** \brief Declarations for the \f$\alpha\f$-QSS method

\author Rowan J Gollan
\version 21-Feb-2006

This method implements the \f$\alpha\f$-QSS method
by Mott.  This method was developed to handle stiff
chemistry integration problems.

**/

class QssStep : public OdeStep {
public:
    /// \brief Normal constructor
    QssStep( const std::string name, int ndim, int max_correctors=4,
	     double qss_eps1=1.1e-5, double c=1.1 );
    
    /// \brief Copy constructor
    QssStep( const QssStep &q );

    /// \brief Default destructor
    virtual ~QssStep();

    /// \brief Clone of the QssStep object
    QssStep* clone();

    // -------- Behaviour of a modified Euler step ------ //

    /// \brief Takes a step using the \f$\alpha\f%-QSS method
    bool advance( OdeSystem &ode, const std::valarray<double> &yin, 
		  std::valarray<double> &yout, double *h );

private:
    int max_correctors_;
    double qss_eps1_;
    double c_;
    double qss_eps2_;
    // Some working arrays for the model
    std::valarray<double> p0_;
    std::valarray<double> q0_;
    std::valarray<double> L0_;
    std::valarray<double> L_p_;
    std::valarray<double> a0_;
    std::valarray<double> y_p_;
    std::valarray<double> y_c_;
    std::valarray<double> p_p_;
    std::valarray<double> p_bar_;
    std::valarray<double> a_bar_;
    std::valarray<double> q_p_;
    std::valarray<double> q_til_;

    // Some private work functions
    void p_on_y( const std::valarray<double> &p, const std::valarray<double> &y,
		 std::valarray<double> &p_y );
    void alpha( double h, const std::valarray<double> &p, std::valarray<double> &a );
    void p_bar( const std::valarray<double> &p0, const std::valarray<double> &p_p, 
		std::valarray<double> &p_bar );
    void q_tilde( const std::valarray<double> &q0, const std::valarray<double> &q_p,
		  const std::valarray<double> &a_bar, std::valarray<double> &q_til );
    bool test_converged( const std::valarray<double> &y_c, const std::valarray<double> &y_p );
    double step_suggest( double h, const std::valarray<double> &y_c,
			 const std::valarray<double> &y_p );

};


/** \brief Declarations for the SAIM method

\author Rowan J Gollan
\version 21-Feb-2006

This class implements the Selected Asymptotic Integrator
developed by Young and Boris.

**/

//------------ Do this! ---------------- //



#endif // ODESTEP_HH
