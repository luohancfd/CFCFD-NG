/** \file segmented_curve_fit.cxx
 * \ingroup nm
 * \brief Definitions of an abtract class for a collection of data curve fits..
 * \author RJG
 * \version 09-Feb-06
 * 
 **/

#include <utility>
#include <vector>
#include <iostream>
#include <sstream>
#include "segmented_curve_fit.hh"
using namespace std;

// ----------------------------------------------------------------------
/// \brief Normal intended constructor for SegmentedCurveFit object.
/// \author RJG
/// \version 10-Feb-06

SegmentedCurveFit::SegmentedCurveFit( const vector<CurveFit*> &segments )
    : CurveFit(), seg_(segments)
{
    double min_val = segments[0]->lower_limit();
    min_index_ = 0;
    double max_val = segments[0]->upper_limit();
    max_index_ = 0;

    for( size_t i = 1; i < segments.size(); ++i) {
	if( segments[i]->lower_limit() < min_val ) {
	    min_val = segments[i]->lower_limit();
	    min_index_ = int(i);
	}
	if ( segments[i]->upper_limit() > max_val ) {
	    max_val = segments[i]->upper_limit();
	    max_index_ = int(i);
	}
    }
    
    lower_limit_ = min_val;
    upper_limit_ = max_val;

    // The initialization above only gave us pointers to segments
    // but we have no control over the persistence of those segments.
    // So here we create new pointers to our own area of memory.

    for( size_t i = 0; i < segments.size(); ++i ) {
	seg_[i] = segments[i]->clone();
    }
    
}


/// \brief Copy constructor for SegmentedCurveFit object
/// \author RJG
/// \version 10-Feb-06

SegmentedCurveFit::SegmentedCurveFit( const SegmentedCurveFit &c )
    : CurveFit( c.lower_limit_, c.upper_limit_ ), seg_( c.seg_ ),
      min_index_( c.min_index_ ), max_index_( c.max_index_ )
{

    for( size_t i = 0; i < c.seg_.size(); ++i ) {
	seg_[i] = c.seg_[i]->clone();
    }

}


/// \brief Default destructor for SegmentedCurveFit object
/// \author RJG
///.\version 10-Feb-06

SegmentedCurveFit::~SegmentedCurveFit()
{
    // Destroy the cloned objects.
    for ( size_t i = 0; i < seg_.size(); ++i ) {
 	delete seg_[i];
    }

}

SegmentedCurveFit* SegmentedCurveFit::clone()
{
    return new SegmentedCurveFit(*this);
}


/// \brief String representation for SegmentedCurveFit object
string SegmentedCurveFit::str() const
{
    ostringstream ost;
    ost << "SegmentedCurveFit( "
	<< "lower_limit=" << lower_limit_
	<< ", upper_limit=" << upper_limit_
	<< "\nwith vector of CurveFits...\n";
    for( vector<CurveFit*>::const_iterator it = seg_.begin();
	 it != seg_.end(); ++it ) {
	ost << **it;
    }
    ost << " )\n";
    return ost.str();
}


double SegmentedCurveFit::eval( double val )
{
 
    if( ! in_range( val ) ) {
	return ( eval_out_of_range( val ) );
    }
    else {
	for( vector<CurveFit*>::iterator it = seg_.begin();
	     it != seg_.end(); ++it ) {
	    if( (*it)->in_range( val ) )
		 return ( (*it)->eval( val ) );
	}
    }
    // We should never reach this point.  But we'll
    // return a value of 0.0.  This can't be handled
    // at a higher level because a value of 0.0
    // might be completely valid.
    // Perhaps this is a candidate to throw and exception.
    return 0.0;
}

double SegmentedCurveFit::eval_out_of_range( double val )
{
    // All we are required to do is defer the work to the segment closest to 
    // the out of range value.
    if( val < lower_limit_ ) {
	return ( seg_[min_index_]->eval_out_of_range( val ) );
    }
    else { // It must be the upper value limit.
	return ( seg_[max_index_]->eval_out_of_range( val ) );
    }

}

// -----------------------------------------------------------------------
// Helper functions for the CurveFit class.
//
// Overload stream output for CurveFit objects
ostream& operator<<( ostream &os, const SegmentedCurveFit &c )
{
    os << c.str();
    return os;
}
