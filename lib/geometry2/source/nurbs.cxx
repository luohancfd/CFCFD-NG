// Author: Rowan J. Gollan
// Date: 10-Sep-2008

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "nurbs.hh"

using namespace std;

//****f*
int find_span(double u, int n, int p, const vector<double> &U)
/** Determines the knot span index.
 *
 *  Input:
 *  u -- find the span in which u lies
 *  n -- there are n+1 basis functions
 *  p -- degree
 *  U -- the knot vector
 *
 *  Return:
 *  i -- the knot span index
 *
 *  This code is an implementation of the C-like pseuocode
 *  given as ALGORITHM 2.1 in Piegl and Tiller (1997). 
 **/
{
    const double TOL = 1.0e-12;

    // Do sanity check first.
    if ( u > (1.0+TOL)*U[n+1] || u < U[0] ) {
	cout << "find_span():\n";
	cout << "u = " << u << " is out of range.\n";
	cout << "U_low = " << U[0] << " U_high= " << U[n+1] << endl;
	cout << "Bailing out!\n";
	exit(1);
    }

    //    if ( u == U[n+1] )
    //	return n;
    // Above is Piegl and Tiller's line, but this wont' work 
    // well for floating point.
    if ( fabs(u - U[n+1]) < TOL )
	return n;

    int low = p;
    int high = n+1;
    int mid = (low + high)/2;
    while ( u < U[mid] || u >= U[mid+1] ) {
	if ( u < U[mid] )
	    high = mid;
	else
	    low = mid;
	mid = (low + high)/2;
    }
    return mid;
}

//****f*
int basis_funs(double u, int i, int p,
	       const vector<double> &U, vector<double> &N)
/** Compute the nonvanishing basis functions.
 *
 * Input:
 * u -- the parameter at which to evaluate basis functions
 * i -- index of basis function
 * p -- the degree
 * U -- the knot vector
 * N -- the vector to be filled with basis functions evaluated
 *      at u for N[0] ... N[p]
 *
 * Return:
 * 0 -- if successfule
 * N -- changed, values of basis functions stored here
 **/
{
    if ( int(N.size()) != p+1 ) {
	cout << "Problem with dimensions of passed in N vector in\n";
	cout << "basis_funs().  A vector of dimension " << p+1 << " is expected.\n";
	cout << "Instead, a vector of dimension " << N.size() << " was given.\n";
	cout << "Bailing out!\n";
	exit(1);
    }
    
    N[0] = 1.0;
    
    // NOTE: This will be inefficient to claim memory on every
    //       call, whether or not its a performance killer remains to be seen.
    //       Perhaps re-implement with passed in working arrays.
    //       Then the Path object could be responsible for managing
    //       the working arrays.
    vector<double> left(p+1, 0.0);
    vector<double> right(p+1, 0.0);
    double saved = 0.0;
    double temp = 0.0;
    for ( int j = 1; j <= p; ++j ) {
	left[j] = u - U[i+1-j];
	right[j] = U[i+j] - u;
	saved = 0.0;
	for( int r = 0; r < j; ++r ) {
	    temp = N[r]/(right[r+1] + left[j-r]);
	    N[r] = saved + right[r+1]*temp;
	    saved = left[j-r]*temp;
	}
	N[j] = saved;
    }
    return 0;
}

//****f*
double one_basis_fun(int p, const vector<double> U,
		     int i, double u)
{
/** Compute the basis function N_i,p(u)
 *
 * Algorithm A2.4 in Piegl and Tiller (1997)
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 * Input:
 * p    -- degree
 * U    -- knot vector
 * i    -- index of basis function
 * u    -- parameter at which to evaluate
 *
 * Return:
 * Nip  -- value of basis function
 *
 **/

    const double ZERO_TOL = 1.0e-15;
    double saved, temp;
    double Uleft, Uright;
    vector<double> N(p+1, 0.0);

    int m = U.size()-1;

    if ( (i == 0 && fabs(u - U[0]) < ZERO_TOL) ||
	 (i == m - p - 1 && fabs(u - U[m]) < ZERO_TOL) ) {
	// Special cases
	return 1.0;
    }

    if ( u < U[i] || u >= U[i+p+1] ) { // local property
	return 0.0;
    }

    for ( int j = 0; j <= p; ++j ) {
	if ( u >= U[i+j] && u < U[i+j+1]) {
	    N[j] = 1.0;
	}
	else {
	    N[j] = 0.0;
	}
    }
	 
    for ( int k = 1; k <= p; ++k ) {
	if ( N[0] == 0.0 )
	    saved = 0.0;
	else
	    saved = ((u - U[i])*N[0])/(U[i+k] - U[i]);
	
	for ( int j = 0; j < p-k+1; ++j ) {
	    Uleft = U[i+j+1];
	    Uright = U[i+j+k+1];
	    
	    if ( N[j+1] == 0.0 ) {
		N[j] = saved;
		saved = 0.0;
	    }
	    else {
		temp = N[j+1]/(Uright-Uleft);
		N[j] = saved+(Uright-u)*temp;
		saved = (u-Uleft)*temp;
	    }
	}
    }

    return N[0];
}

//****f*
Vector3 nurbs_curve_point(double u, int p,
			  const vector<double> &U,
			  const vector<Mapped_point> &Pw)
/** Compute a point on the rational B-spline curve.
 *
 *  Input:
 *  u -- the parameter at which to evaluate
 *  p -- the degree
 *  U -- the knot vector
 *  Pw -- the vector of mapped control points (weights included)
 *
 *  Return:
 *  C -- point on curve (Vector3)
 **/
{
    int m = U.size() - 1;
    int n = m - p - 1;

    int span = find_span(u, n, p, U);
    vector<double> N(p+1, 0.0);
    basis_funs(u, span, p, U, N);

    Mapped_point Cw;
    Cw.wx = 0.0; Cw.wy = 0.0; Cw.wz = 0.0; Cw.w = 0.0;

    for ( int j = 0; j <= p; ++j ) {
	// This is the evaluation in 4D space.
	Cw = Cw + N[j]*Pw[span-p+j];
    }

    return Cw.to_Vector3();
}

//****f*
Vector3 nurbs_surface_point(double u, int p,
			    const vector<double> &U,
			    double v, int q,
			    const vector<double> &V,
			    const vector<vector<Mapped_point> >&Pw)
/** Compute a point on the rational B-spline surface.
 *
 *  Input:
 *  u -- u parameter at which to evaluate
 *  p -- the degree of u curves
 *  U -- the knot vector for u curves
 *  v -- v parameter at which to evaluate
 *  q -- the degree of v curves
 *  V -- the knot vector for v curves
 *  Pw -- the net of mapped control points (weights included) 
 *
 *  Return:
 *  S -- point on surface (Vector3)
 **/
{
    int m = U.size() - 1;
    int n = m - p - 1;
    int uspan = find_span(u, n, p, U);
    vector<double> Nu(p+1, 0.0);
    basis_funs(u, uspan, p, U, Nu);

    m = V.size() - 1;
    n = m - q - 1;
    int vspan = find_span(v, n, q, V);
    vector<double> Nv(q+1, 0.0);
    basis_funs(v, vspan, q, V, Nv);

    vector<Mapped_point> temp(q+1);

    for ( int l = 0; l <= q; ++l ) {
	for ( int k = 0; k <= p; ++k ) {
	    temp[l] = temp[l] + Nu[k]*Pw[uspan-p+k][vspan-q+l];
	}
    }

    Mapped_point Sw;
    
    for ( int l = 0; l <= q; ++l )
	Sw = Sw + Nv[l]*temp[l];

    return Sw.to_Vector3();
}


Mapped_point::
Mapped_point()
    : wx(0.0), wy(0.0), wz(0.0), w(0.0) {}

Mapped_point::
Mapped_point(double wx, double wy, double wz, double w)
    : wx(wx), wy(wy), wz(wz), w(w) {}

Mapped_point::
Mapped_point(const Mapped_point &mp)
    : wx(mp.wx), wy(mp.wy), wz(mp.wz), w(mp.w) {}

Mapped_point::
~Mapped_point() {}

Mapped_point&
Mapped_point::
operator=(const Mapped_point &mp)
{
    if( this == &mp ) // avoid aliasing
	return (*this);
    
    wx = mp.wx;
    wy = mp.wy;
    wz = mp.wz;
    w = mp.w;

    return (*this);
}

Vector3
Mapped_point::
to_Vector3()
{
    Vector3 v(wx, wy, wz);
    if ( w != 0.0 ) {
	v = v/w;
    }
    // else infinite control point
    // do nothing.
    
    return v;
}

Mapped_point operator+(const Mapped_point &m1, const Mapped_point &m2)
{
    return Mapped_point(m1.wx + m2.wx,
			m1.wy + m2.wy,
			m1.wz + m2.wz,
			m1.w + m2.w);
}

Mapped_point operator-(const Mapped_point &m1, const Mapped_point &m2)
{
    return Mapped_point(m1.wx - m2.wx,
			m1.wy - m2.wy,
			m1.wz - m2.wz,
			m1.w - m2.w);
}

Mapped_point operator*(double scalar, const Mapped_point &mp)
{
    return Mapped_point(scalar*mp.wx,
			scalar*mp.wy,
			scalar*mp.wz,
			scalar*mp.w);
}

Mapped_point operator*(const Mapped_point &mp, double scalar)
{
    return Mapped_point(scalar*mp.wx,
			scalar*mp.wy,
			scalar*mp.wz,
			scalar*mp.w);
}

double abs(const Mapped_point &m)
{
    return sqrt(m.wx*m.wx + m.wy*m.wy + m.wz*m.wz + m.w*m.w);
}
