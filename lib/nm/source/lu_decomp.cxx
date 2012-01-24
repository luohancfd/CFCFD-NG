// Author: Rowan J. Gollan
// Date: 02-Mar-2010
// Place: Poquoson, Virginia, USA
// 
// See header file for notes.
//

#include <cmath>

#include "lu_decomp.hh"

using namespace std;

LUdcmp::
LUdcmp(const Valmatrix &a)
    : n_(a.nrows()), lu_(a), indx_(n_, 0)
{
    const double TINY = 1.0e-40;
    size_t i, j, k;
    size_t imax = 0;
    double big, temp;
    vector<double> vv(n_, 0.0);

    d_ = 1.0;
    
    for ( i = 0; i < n_; ++i ) {
	big = 0.0;
	for ( j = 0; j < n_; ++j ) {
	    if ( (temp=fabs(lu_.get(i,j))) > big )
		big = temp;
	}
	if ( big == 0.0 )
	    throw("Singular matrix in LUdcmp");

	vv[i] = 1.0/big;
    }

    for ( k = 0; k < n_; ++k ) {
	big = 0.0;
	for ( i = k; i < n_; ++i ) {
	    temp = vv[i] * fabs(lu_.get(i,k));
	    if ( temp > big ) {
		big = temp;
		imax = i;
	    }
	}
	if ( k != imax ) {
	    for ( j = 0; j < n_; ++j ) {
		temp = lu_.get(imax,j);
		lu_.set(imax,j, lu_.get(k,j));
		lu_.set(k,j,temp);
	    }
	    d_ = -d_;
	    vv[imax] = vv[k];
	}
	indx_[k] = imax;
	if ( lu_.get(k,k) == 0.0 )
	    lu_.set(k,k,TINY);

	for ( i = k+1; i < n_; ++i ) {
	    lu_.set(i,k, lu_.get(i,k)/lu_.get(k,k));
	    temp = lu_.get(i,k);
	    for ( j = k+1; j < n_; ++j ) {
		lu_.set(i,j, lu_.get(i,j) - temp*lu_.get(k,j));
	    }
	}
    }
}

void
LUdcmp::
solve(const vector<double> &b, vector<double> &x)
{
    size_t i, ip, j;
    size_t ii = 0;
    double sum;

    if ( (b.size() != n_) || (x.size() != n_) )
	throw("LUdcmp::solve bad sizes");

    for ( i = 0; i < n_; ++i )
	x[i] = b[i];

    for ( i = 0; i < n_; ++i ) {
	ip = indx_[i];
	sum = x[ip];
	x[ip] = x[i];
	if ( ii != 0 )
	    for ( j = ii-1; j < i; ++j ) sum -= lu_.get(i,j)*x[j];
	else if ( sum != 0.0 )
	    ii = i + 1;
	x[i] = sum;
    }

    // Need an integer here to do back-substitution
    int I;
    for ( I = n_ - 1; I >= 0; --I ) {
	sum = x[I];
	for ( j = I+1; j < n_; ++j ) sum -= lu_.get(I,j)*x[j];
	x[I] = sum/lu_.get(I,I);
	    
    }
}

