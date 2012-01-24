/** \file curve_fit.hh
 * \ingroup nm
 * \brief An abstract class for data curve fits.
 * \author RJG
 * \version 09-Feb-06
 * 
 **/

#ifndef CURVE_FIT_HH
#define CURVE_FIT_HH

#include <utility>
using namespace std;

// -----------------------------------------------------------------------------
/// \brief An abstract class defining the public face of a curve fit.
/// \author RJG
/// \version 09-Feb-06
///
/// In defining the abstract base class, we ask some questions:
/// Q. What are the operations common to all curve fits?
/// A. All curve fits need a publicly accessible eval() function.
///    Furthermore, they need to be able to report that an input value
///    is in the range of validity for the curve fit.
/// Q. Follow on question: what minimal set of data is required?
/// A. We may not know the form of the data required to represent the curve fit
///    but we do know that it will have a range of validity represented by a 
///    pair of doubles.
///
 
class CurveFit {
public:
    // Normal contructor
    CurveFit( double low_val, double high_val );
    CurveFit();

    // Copy constructor
    CurveFit( const CurveFit &c );

    // Destructor
    virtual ~CurveFit();

    virtual CurveFit* clone() = 0;

    // The publicly exposed eval function.  This is the real
    // workhorse of the curve fit and is the essential behaviour
    // we hope to encapsulate with this class.
    virtual double eval( double val ) = 0;

    /// \brief A pure virtual eval_out_of_range function.
    /// \author RJG
    /// \version 09-Feb-06
    /// This is a pure virtual function indicating that derived
    /// classes of CurveFit must provide their own implementation.
    virtual double eval_out_of_range( double val ) = 0;

    // The publicly exposed in_range function returns a Boolean
    // letting the caller know if their value is in the range of
    // validity for the curve fit.
    bool in_range( double val );

    // Let the outside world know the range if requested.
    double lower_limit() { return lower_limit_; }
    double upper_limit() { return upper_limit_; }

    /// \brief Returns a string representation of the CurveFit.
    virtual string str() const;

protected:
    // Data available in the derived classes BUT not to the caller
    double lower_limit_;
    double upper_limit_;

};

// Helper functions.

/// Output a string representation of the generic curve fit.
ostream& operator<<( ostream &os, const CurveFit &c );

#endif
