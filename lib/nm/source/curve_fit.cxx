/** \file curve_fit.cxx
 * \ingroup nm
 * \brief Definitions of an abtract class for data curve fits..
 * \author RJG
 * \version 09-Feb-06
 * 
 **/

#include <utility>
#include <iostream>
#include <sstream>
#include "curve_fit.hh"
using namespace std;

// ----------------------------------------------------------------------
/// \brief Normal intended constructor for CurveFit object.
/// \author RJG
/// \version 09-Feb-06

CurveFit::CurveFit( double low_val, double high_val )
    : lower_limit_( low_val ), upper_limit_( high_val ) {}

CurveFit::CurveFit() {}

/// \brief Copy constructor for CurveFit object
/// \author RJG
/// \version 09-Feb-06

CurveFit::CurveFit( const CurveFit &c ) {}

/// \brief Default destructor for CurveFit object
/// \author RJG
///.\version 09-Feb-06

CurveFit::~CurveFit() {}


/// \brief String representation for CurveFit object
string CurveFit::str() const
{
    ostringstream ost;
    ost << "CurveFit( "
	<< "lower_limit=" << lower_limit_
	<< ", high_val=" << upper_limit_
	<< " )\n";
    return ost.str();
}

/// \brief The in_range function common to all CurveFits.
/// \author RJG
/// \version 09-Feb-06

bool CurveFit::in_range( double val )
{
    return ( val >= lower_limit_ && val <= upper_limit_ );
}
    
// -----------------------------------------------------------------------
// Helper functions for the CurveFit class.
//
// Overload stream output for CurveFit objects
ostream& operator<<( ostream &os, const CurveFit &c )
{
    os << c.str();
    return os;
}
