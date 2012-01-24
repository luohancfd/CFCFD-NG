/** \file segmented_curve_fit.hh
 * \ingroup nm
 * \brief An abstract class for a collection of data curve fits.
 * \author RJG
 * \version 09-Feb-06
 * 
 **/

#ifndef SEG_CURVE_FIT_HH
#define SEG_CURVE_FIT_HH

#include <utility>
#include <vector>
#include <iostream>
#include <sstream>
#include "curve_fit.hh"

using namespace std;

// -----------------------------------------------------------------------------
/// \brief An abstract class defining the public face of a segmented curve fit.
/// \author RJG
/// \version 09-Feb-06
///
/// This is a derived class of CurveFit but it is still abstract in that 
/// it doesn't represent any particular CurveFit.
///
 
class SegmentedCurveFit : public CurveFit {
public:
    // Normal contructor
    SegmentedCurveFit( const vector<CurveFit*> &segments );

    // Copy constructor
    SegmentedCurveFit( const SegmentedCurveFit &c );

    // Destructor
    virtual ~SegmentedCurveFit();

    virtual SegmentedCurveFit* clone();

    // The publicly exposed eval function.
    // This function is the same in all derived classes
    // because the segmented curve fit always behaves 
    // in the same way. The underlying collection of
    // of segments may behave as their subclass dictates.
    double eval( double val );

     /// \brief A pure virtual eval_out_of_range function.
    /// \author RJG
    /// \version 09-Feb-06
    /// This is a pure virtual function indicating that derived
    /// classes of SegmentedCurveFit must provide their own implementation.
    virtual double eval_out_of_range( double val );

    /// \brief Returns a string representation of the CurveFit.
    virtual string str() const;

protected:
    // The extra data in a SegmentedCurveFit is the collection
    // of CurevFit segments.
    vector<CurveFit*> seg_;
    int min_index_;
    int max_index_;

};

/// Output a string representation of the generic segmented curve fit.
ostream& operator<<( ostream &os, const SegmentedCurveFit &c );

#endif


