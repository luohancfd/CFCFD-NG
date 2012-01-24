// Author: Rowan J. Gollan
// Date: 04-Nov-2008
// Place: Hampton, Virginia, USA

#include <vector>

#include "segmented-functor.hh"

using namespace std;

Segmented_functor::
Segmented_functor(vector<Univariate_functor*> &f, vector<double> &breaks)
    : breaks_(breaks)
{
    for ( size_t i = 0; i < f.size(); ++i ) {
	f_.push_back(f[i]->clone());
    }
}

Segmented_functor::
Segmented_functor(const Segmented_functor &s)
    : breaks_(s.breaks_)
{
    for ( size_t i = 0; i < s.f_.size(); ++i ) {
	f_.push_back(s.f_[i]->clone());
    }
}

Segmented_functor::
~Segmented_functor()
{
    for ( size_t i = 0; i < f_.size(); ++i ) {
	delete f_[i];
    }
}

Segmented_functor*
Segmented_functor::
clone() const
{
    return new Segmented_functor(*this);
}


double
Segmented_functor::
operator()(double x)
{
    if ( x < breaks_[0] ) {
	// Let the first functor take care of the result
	return (*f_[0])(x);
    }
    
    for ( size_t i = 0; i < f_.size(); ++i ) {
	// Loop through looking for appropriate functor
	if ( x >= breaks_[i] && x < breaks_[i+1] )
	    return (*f_[i])(x);
    }
    
    // If we reached this point, just let the last functor
    // take care of the evaluation
    return (*f_.back())(x);
}
