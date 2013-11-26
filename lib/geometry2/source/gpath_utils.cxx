/** \file gpath_utils.cxx
 *  \ingroup libgeom2
 *  \brief Definitions for the C++ geometric-path utitlities.
 *  \author RJG
 *  \version 07-Feb-2008 -- intial code
 *
 */

#include <iostream>
#include <vector>
#include <cmath>

#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "../../nm/source/fobject.hh"
#include "../../nm/source/nelmin.hh"

#include "geom.hh"
#include "gpath.hh"
#include "gpath_utils.hh"

using namespace std;

#define A(i) ((i)*6)
#define B(i) ((i)*6+1)
#define C(i) ((i)*6+2)
#define D(i) ((i)*6+3)
#define E(i) ((i)*6+4)
#define F(i) ((i)*6+5)

//****f*
XSpline create_quintic_spline(vector<Vector3> &P,
			      vector<double> &Pdash,
			      double d2a, double d2b)
{
/** Create a quintic spline.
 *
 * This function returns a quintic spline of type XSpline
 * interpolating a set of points and derivatives at those points.
 *
 * Inputs:
 * P     -- a vector of Vector3s: these are the points for interpolation
 * Pdash -- a vector of doubles: these are the derivatives (dydx)
 *          for interpolation.
 * d2a   -- the second derivative (curvature) at the left end of the spline
 * d2b   -- the second derivative (curvature) at the right end of the spline
 **/

    if( P.size() != Pdash.size() ) {
	cout << "create_quintic_spline(): Sizes of vector of points and derivatives\n";
	cout << "                         do not match.\n";
	cout << "Bailing Out!\n";
	exit(1);
    }

    int offset;
    size_t N = P.size() - 1;
    Valmatrix M(6*N, 6*N, 0.0);
    vector<double> U(6*N);
    
    // Ensure U is filled with zeros.
    for( size_t i = 0; i < U.size(); ++i ) U[i] = 0.0;


    // Populate M and U representing the linear
    // set of equations.

    // For all segments of the spline we know the
    // points and derivatives.  So let's specify this
    // in the equations.
    for( size_t i = 0; i < N; ++i ) {
	offset = i*4;
	// S_i(x_i) = P[i].y
	M.set(offset, A(i), 1.0);
	M.set(offset, B(i), P[i].x);
	M.set(offset, C(i), P[i].x*P[i].x);
	M.set(offset, D(i), pow(P[i].x, 3));
	M.set(offset, E(i), pow(P[i].x, 4));
	M.set(offset, F(i), pow(P[i].x, 5));
	U[offset] = P[i].y;
	// S_i'(x_i) = Pdash[i]
	M.set(offset+1, B(i), 1.0);
	M.set(offset+1, C(i), 2.0*P[i].x);
	M.set(offset+1, D(i), 3.0*P[i].x*P[i].x);
	M.set(offset+1, E(i), 4.0*pow(P[i].x, 3));
	M.set(offset+1, F(i), 5.0*pow(P[i].x, 4));
	U[offset+1] = Pdash[i];
	// S_i(x_(i+1)) = P[i+1].y
	M.set(offset+2, A(i), 1.0);
	M.set(offset+2, B(i), P[i+1].x);
	M.set(offset+2, C(i), P[i+1].x*P[i+1].x);
	M.set(offset+2, D(i), pow(P[i+1].x, 3));
	M.set(offset+2, E(i), pow(P[i+1].x, 4));
	M.set(offset+2, F(i), pow(P[i+1].x, 5));
	U[offset+2] = P[i+1].y;
	// S_i'(x_i+1) = Pdash[i+1]
	M.set(offset+3, B(i), 1.0);
	M.set(offset+3, C(i), 2.0*P[i+1].x);
	M.set(offset+3, D(i), 3.0*P[i+1].x*P[i+1].x);
	M.set(offset+3, E(i), 4.0*pow(P[i+1].x, 3));
	M.set(offset+3, F(i), 5.0*pow(P[i+1].x, 4));
	U[offset+3] = Pdash[i+1];
    }
    
    // Now the equations related to the continuity
    // of second and third derivatives at the knot points.
    for( size_t i = 0; i < (N-1); ++i ) {
	offset = 4*N + i*2;
	// S_i''(x_(i+1)) = S_(i+1)''(x_(i+1))
	M.set(offset, C(i), 2.0);
	M.set(offset, D(i), 6.0*P[i+1].x);
	M.set(offset, E(i), 12.0*P[i+1].x*P[i+1].x);
	M.set(offset, F(i), 20.0*pow(P[i+1].x, 3));
	M.set(offset, C(i+1), -2.0);
	M.set(offset, D(i+1), -6.0*P[i+1].x);
	M.set(offset, E(i+1), -12.0*P[i+1].x*P[i+1].x);
	M.set(offset, F(i+1), -20.0*pow(P[i+1].x, 3));
	U[offset] = 0.0;
	// S_i'''(x_(i+1)) = S_(i+1)'''(x_(i+1))
	M.set(offset+1, D(i), 6.0);
	M.set(offset+1, E(i), 24.0*P[i+1].x);
	M.set(offset+1, F(i), 60.0*P[i+1].x*P[i+1].x);
	M.set(offset+1, D(i+1), -6.0);
	M.set(offset+1, E(i+1), -24.0*P[i+1].x);
	M.set(offset+1, F(i+1), -60.0*P[i+1].x*P[i+1].x);
	U[offset+1] = 0.0;
    }

    // Now the boundaries.
    // Left end:
    M.set(6*N-2, C(0), 2.0);
    M.set(6*N-2, D(0), 6.0*P[0].x);
    M.set(6*N-2, E(0), 12.0*P[0].x*P[0].x);
    M.set(6*N-2, F(0), 20.0*pow(P[0].x, 3));
    U[6*N-2] = d2a;
    // Right end:
    M.set(6*N-1, C(N-1), 2.0);
    M.set(6*N-1, D(N-1), 6.0*P[N].x);
    M.set(6*N-1, E(N-1), 12.0*P[N].x*P[N].x);
    M.set(6*N-1, F(N-1), 20.0*pow(P[N].x, 3));
    U[6*N-1] = d2b;

    // END setting up equations.

    vector<double> x(6*N);
    gaussian_elimination(M, x, U, false);

    vector<XPoly> polys;
    for( size_t i = 0; i < N; ++i ) {
	vector<double> tmp;
	tmp.push_back(x[A(i)]);
	tmp.push_back(x[B(i)]);
	tmp.push_back(x[C(i)]);
	tmp.push_back(x[D(i)]);
	tmp.push_back(x[E(i)]);
	tmp.push_back(x[F(i)]);
	polys.push_back(XPoly(P[i].x, P[i+1].x, tmp));
    }
    return XSpline(P[0].x, P[N].x, polys);
}

#undef A
#undef B
#undef C
#undef D
#undef E
#undef F

#define A(i) ((i)*4)
#define B(i) ((i)*4+1)
#define C(i) ((i)*4+2)
#define D(i) ((i)*4+3)

//****f*
XSpline create_cubic_spline(vector<Vector3> &P,
			      vector<double> &Pdash,
			      double d2a, double d2b)
{
/** Create a cubic spline.
 *
 * This function returns a cubic spline of type XSpline
 * interpolating a set of points and derivatives at those points.
 *
 * This interpolation is a "best-fit" type interpolation: the
 * equation set is overdetermined and so a linear least squares
 * solution is found.
 *
 * Inputs:
 * P     -- a vector of Vector3s: these are the points for interpolation
 * Pdash -- a vector of doubles: these are the derivatives (dydx)
 *          for interpolation.
 * d2a   -- the second derivative (curvature) at the left end of the spline
 * d2b   -- the second derivative (curvature) at the right end of the spline
 **/
    if( P.size() != Pdash.size() ) {
	cout << "create_cubic_spline(): Sizes of vector of points and derivatives\n";
	cout << "                       do not match.\n";
	cout << "Bailing Out!\n";
	exit(1);
    }

    double del_x;
    int offset;
    size_t N = P.size() - 1;
    Valmatrix M(5*N+1, 4*N, 0.0);
    vector<double> U(5*N+1);

    // Ensure U is filled with zeros.
    for( size_t i = 0; i < U.size(); ++i ) U[i] = 0.0;
    // Populate M and U representing the linear
    // set of equations.

    // For all segments of the spline we know the
    // points and derivatives.  So let's specify this
    // in the equations.
    for( size_t i = 0; i < N; ++i ) {
	del_x = P[i+1].x - P[i].x;
	offset = i*4;
	// S_i(x_i) = P[i].y
	M.set(offset, A(i), 1.0);
	U[offset] = P[i].y;
	// S_i'(x_i) = Pdash[i]
	M.set(offset+1, B(i), 1.0);
	U[offset+1] = Pdash[i];
	// S_i(x_i+1) = P[i+1].y
	M.set(offset+2, A(i), 1.0);
	M.set(offset+2, B(i), del_x);
	M.set(offset+2, C(i), del_x*del_x);
	M.set(offset+2, D(i), pow(del_x, 3));
	U[offset+2] = P[i+1].y;
	// S_i'(x_i+1) = Pdash[i+1]
	M.set(offset+3, B(i), 1.0);
	M.set(offset+3, C(i), 2.0*del_x);
	M.set(offset+3, D(i), 3.0*del_x*del_x);
	U[offset+3] = Pdash[i+1];
    }
    // Now the equations related to the continuity
    // of the second derivatives at the knot points
    for( size_t i = 0; i < (N-1); ++i ) {
	del_x = P[i+1].x - P[i].x;
	offset = 4*N + i;
	// S_i''(x_i+1) = S_i+1''(x_i+1)
	M.set(offset, C(i), 2.0);
	M.set(offset, D(i), 6.0*del_x);
	M.set(offset, C(i+1), -2.0);
	U[offset] = 0.0;
    }

    // Finally the curvature at the boundaries...
    // Left end
    M.set(5*N-1, C(0), 2.0);
    U[5*N-1] = d2a;
    // Right end
    del_x = P[N].x - P[N-1].x;
    M.set(5*N, C(N-1), 2.0);
    M.set(5*N, D(N-1), 6.0*del_x);
    U[5*N] = d2b;
    // END of setting up equations.
    
    vector<double> x(4*N);
    least_squares_solve(M, x, U);
    
    vector<XPoly> polys;
    for( size_t i = 0; i < N; ++i ) {
	vector<double> tmp;
	tmp.push_back(x[A(i)]);
	tmp.push_back(x[B(i)]);
	tmp.push_back(x[C(i)]);
	tmp.push_back(x[D(i)]);
	polys.push_back(XPoly(P[i].x, P[i+1].x, tmp));
    }
    return XSpline(P[0].x, P[N].x, polys);
}

#undef A
#undef B
#undef C
#undef D

XPoly create_line(Vector3 &P, double m, double x0, double x1)
{
    vector<double> B;
    B.resize(2);

    // Set gradient
    B[1] = m;
    // Compute B[0]
    B[0] = P.y - m*P.x;
    
    return XPoly(x0, x1, B);
}

int fac(int n)
{
    if( n == 0 )
	return 1;
    return n*fac(n-1);
}

double choose_func(int n, int i)
{
    return double( fac(n) ) / double ( fac(i) * fac(n-i) );
}

double Bernstein(int n, int i, double t)
{
    return choose_func(n, i)*pow(t, i)*pow(1.0 - t, (n-i));
}

//****f*
XBezier create_best_fit_XBezier(vector<Vector3> &P, int n, double lslope, double rslope,
				double lcurv, double rcurv)
{
/** Create a best fit XBezier to a set of data points.
 *
 * This function finds the best fit Bezier of specified order
 * to a set of data points (x,y).  It also meets the constraints
 * of specified slope and curvature (approximate) at the end points.
 *
 * Being an XBezier, the control points will be spaced
 * equidistantly in the x-direction.  This does limit how well the
 * the fit is.
 *
 * Inputs:
 * P     -- a vector of Vector3s: date points for approximation
 * n     -- order of Bezier curve
 * lslope -- slope at left end of curve (x = x0)
 * rslope -- slope at right end of curve (x = x1)
 * lcurv  -- (approximate) curvature at left end of curve
 * rcurv  -- (approximate) curvature at right end of curve
 *
 * By "approximate" curvature, it is meant that curvature
 * is treated simply as y'' rather than:
 * Kappa = y''/(1 + (y')^2)^(2/3)
 *
 * For lcurv = rcurv = 0.0, this approximate treatment is exact.
 *
 **/
    int nc = n + 1; // no. control npoints
    int nu = nc - 6; // no. unknowns
    int Np = P.size() - 2; // no. points to fit
                           // we do not need to "fit" the end points,
                           // hence minus 2
    vector<Vector3> B;
    B.resize(nc);
    B[0] = P[0];
    B.back() = P.back();
    
    double x0 = P[0].x;
    double x1 = P.back().x;
    
    // Set the aspects of the B-coefficients 
    // which are set by specified slope and curvature

    // B[1] set by lslope
    B[1].y = (x1 - x0)*lslope/n + B[0].y;
    // B[2] set by lcurv
    B[2].y = (x1 - x0)*lcurv/(n*(n-1)) + 2.0*B[1].y - B[0].y;
    // B[-2] set by rslope
    B[B.size()-2].y = B.back().y - (x1 - x0)*rslope/n;
    // B[-3] set by rcurv
    B[B.size()-3].y = (x1 - x0)*rcurv/(n*(n-1)) - B.back().y + 2.0*B[B.size()-2].y;

    // The XBezier requires the x locations of control points
    // to be placed regularly.
    double dx = (x1 - x0)/(nc-1);
    for( int i = 0; i < nc; ++i ) B[i].x = x0 + i*dx;

    Valmatrix A(Np, nu, 0.0);
    vector<double> b(Np);

    for( size_t i = 0; i < b.size(); ++i ) b[i] = 0.0;

    double x, t, RHS;

    // Equations in y...
    if( nc >= 6 ) {
	for( size_t i = 1; i < P.size()-1; ++i ) {
	    x = P[i].x;
	    t = (x - x0)/(x1 - x0);
	    RHS = P[i].y - (B[0].y*Bernstein(n, 0, t) + B[1].y*Bernstein(n, 1, t) +
			    B[2].y*Bernstein(n, 2, t) + B[B.size()-3].y*Bernstein(n, n-2, t) +
			    B[B.size()-2].y*Bernstein(n, n-1, t) + B.back().y*Bernstein(n, n, t));
	    for( size_t j = 3; j < size_t(n-2); ++j ) {
		A.set(i-1, j-3, Bernstein(n, j, t));
	    }
	    b[i-1] = RHS;
	}
    }
    
    vector<double> u(nu);
    least_squares_solve(A, u, b);

    if( nc >= 6 ) {
	for( size_t j = 3; j < size_t(n-2); ++j ) {
	    B[j].y = u[j-3];
	}
    }
    
    return XBezier(B[0].x, B.back().x, B);
}

//****f*
Bezier create_optimised_Bezier(vector<Vector3> &P, int n,
			       double dydxA, double dydxB, double z_pos)
{
/** Create a best fit Bezier to a set of data points.
 *
 * This function finds the best fit Bezier of specified order
 * to a set of data points (x,y).  It also meets the constraint
 * of specified slope at the end points.
 *
 * Inputs:
 * P      -- a vector of Vector3s: data points for approximation
 * n      -- order of Bezier curve
 * dydxA -- slope at start of curve (t = 0.0)
 * dydxB -- slope at end of curve (t = 0.0)
 *
 **/
    MultivariateFunction *bfb = new Best_fit_Bezier(P, n, dydxA, dydxB, z_pos);
    
    double f_min = 0;
    int n_fe = 0;
    int n_restart = 0;
    
    // Setup initial guess. A line between first and last point.
    int nc = n+1;
    vector<double> x;
    x.resize((nc - 2) + (nc - 4)); // unknowns (in x, and in y)
    double m = (P.back().y - P[0].y)/(P.back().x - P[0].x);
    double dx = (P.back().x - P[0].x)/(nc - 1);
    x[0] = P[0].x + dx;
    for ( int i = 2; i < nc-2; ++i ) {
	x[2*(i-1)-1] = P[0].x + i*dx;
	x[2*(i-1)] = m*(x[2*(i-1)-1] - P[0].x) + P[0].y;
    }
    x.back() = P.back().x - dx;

    // Call optimiser
    minimize(bfb, x, &f_min, &n_fe, &n_restart);

    // Setup Bezier.
    vector<Vector3> B;
    B.resize(nc);
    B[0] = P[0];
    B.back() = P.back();
    B[1].x = x[0];
    B[1].y = dydxA*(x[0] - B[0].x) + B[0].y;

    B[B.size()-2].x = x.back();
    B[B.size()-2].y = dydxB*(x.back() - B.back().x) + B.back().y;

    for ( size_t i = 2; i < B.size()-2; ++i ) {
	B[i].x = x[2*(i-1)-1];
	B[i].y = x[2*(i-1)];
    }

    // Flush small non-zero values to zero.
    const double small_val = 1.0e-16;
    for ( size_t i = 0; i < B.size(); ++i ) {
	if ( fabs(B[i].x) <= small_val )
	    B[i].x = 0.0;
	if ( fabs(B[i].y) <= small_val )
	    B[i].y = 0.0;
    }

    // Set all B[].z to z_pos
    for ( size_t i = 0; i < B.size(); ++i ) {
	B[i].z = z_pos;
    }

    delete bfb;
    return Bezier(B);
}

//****f*
Bezier create_optimised_Bezier_YZ(vector<Vector3> &P, int n,
				  double dydzA, double dydzB, double x_pos)
{
/** Create a best fit Bezier to a set of data points.
 *
 * This function finds the best fit Bezier of specified order
 * to a set of data points (z,y).  It also meets the constraint
 * of specified slope at the end points.
 *
 * Inputs:
 * P      -- a vector of Vector3s: data points for approximation
 * n      -- order of Bezier curve
 * dydzA -- slope at start of curve (t = 0.0)
 * dydzB -- slope at end of curve (t = 0.0)
 *
 **/
    MultivariateFunction *bfb = new Best_fit_Bezier_YZ(P, n, dydzA, dydzB, x_pos);
    
    double f_min = 0;
    int n_fe = 0;
    int n_restart = 0;
    
    // Setup initial guess. A line between first and last point.
    int nc = n+1;
    vector<double> z;
    z.resize((nc - 2) + (nc - 4)); // unknowns (in z, and in y)
    double m = (P.back().y - P[0].y)/(P.back().z - P[0].z);
    double dz = (P.back().z - P[0].z)/(nc - 1);
    z[0] = P[0].z + dz;
    for ( int i = 2; i < nc-2; ++i ) {
	z[2*(i-1)-1] = P[0].z + i*dz;
	z[2*(i-1)] = m*(z[2*(i-1)-1] - P[0].z) + P[0].y;
    }
    z.back() = P.back().z - dz;

    // Call optimiser
    minimize(bfb, z, &f_min, &n_fe, &n_restart);

    // Setup Bezier.
    vector<Vector3> B;
    B.resize(nc);
    B[0] = P[0];
    B.back() = P.back();
    B[1].z = z[0];
    B[1].y = dydzA*(z[0] - B[0].z) + B[0].y;

    B[B.size()-2].z = z.back();
    B[B.size()-2].y = dydzB*(z.back() - B.back().z) + B.back().y;

    for ( size_t i = 2; i < B.size()-2; ++i ) {
	B[i].z = z[2*(i-1)-1];
	B[i].y = z[2*(i-1)];
    }

    // Flush small non-zero values to zero.
    const double small_val = 1.0e-16;
    for ( size_t i = 0; i < B.size(); ++i ) {
	if ( fabs(B[i].z) <= small_val )
	    B[i].z = 0.0;
	if ( fabs(B[i].y) <= small_val )
	    B[i].y = 0.0;
    }

    // Set all B[].x to x_pos
    for ( size_t i = 0; i < B.size(); ++i ) {
	B[i].x = x_pos;
    }

    delete bfb;
    return Bezier(B);
}


//****f*
Bezier create_optimised_Bezier3D(vector<Vector3> &P, int n,
				 double dydxA, double dydxB,
				 double dzdxA, double dzdxB)
{
/** Create a best fit Bezier (3D) to a set of data points.
 *
 * This function finds the best fit Bezier of specified order
 * to a set of data points (x,y,z).  It also meets the constraint
 * of specified slope at the end points.
 *
 * Inputs:
 * P      -- a vector of Vector3s: date points for approximation
 * n      -- order of Bezier curve
 * dydxA -- dydx slope at start of curve (t = 0.0)
 * dydxB -- dydx slope at end of curve (t = 0.0)
 * dzdxA -- dzdx slope at start of curve
 * dzdxB -- dzdx slope at end of curve
 *
 **/
    MultivariateFunction *bfb = new Best_fit_Bezier3D(P, n, dydxA, dydxB,
						      dzdxA, dzdxB);
    double f_min = 0;
    int n_fe = 0;
    int n_restart = 0;
    
    // Setup initial guess. A line between first and last point.
    int nc = n+1;
    vector<double> x;
    x.resize((nc - 2) + (nc - 4) + (nc - 4)); // unknowns (in x, and in y)
    double m1 = (P.back().y - P[0].y)/(P.back().x - P[0].x);
    double m2 = (P.back().z - P[0].z)/(P.back().x - P[0].x);
    double dx = (P.back().x - P[0].x)/(nc - 1);
    x[0] = P[0].x + dx;
    for ( int i = 2; i < nc-2; ++i ) {
	x[3*(i-1)-2] = P[0].x + i*dx;
	x[3*(i-1)-1] = m1*(x[3*(i-1)-2] - P[0].x) + P[0].y;
	x[3*(i-1)] = m2*(x[3*(i-1)-2] - P[0].x) + P[0].z;
    }
    x.back() = P.back().x - dx;

    // Call optimiser
    minimize(bfb, x, &f_min, &n_fe, &n_restart);

    // Setup Bezier.
    vector<Vector3> B;
    B.resize(nc);
    B[0] = P[0];
    B.back() = P.back();
    B[1].x = x[0];
    B[1].y = dydxA*(x[0] - B[0].x) + B[0].y;
    B[1].z = dzdxA*(x[0] - B[0].x) + B[0].z;

    B[B.size()-2].x = x.back();
    B[B.size()-2].y = dydxB*(x.back() - B.back().x) + B.back().y;
    B[B.size()-2].z = dzdxB*(x.back() - B.back().x) + B.back().z;

    for ( size_t i = 2; i < B.size()-2; ++i ) {
	B[i].x = x[3*(i-1)-2];
	B[i].y = x[3*(i-1)-1];
	B[i].z = x[3*(i-1)];
    }

    // Flush small non-zero values to zero.
    const double small_val = 1.0e-16;
    for ( size_t i = 0; i < B.size(); ++i ) {
	if ( fabs(B[i].x) <= small_val )
	    B[i].x = 0.0;
	if ( fabs(B[i].y) <= small_val )
	    B[i].y = 0.0;
	if ( fabs(B[i].z) <= small_val )
	    B[i].z = 0.0;
    }

    delete bfb;
    return Bezier(B);
}

//****f*
XPoly create_best_fit_poly(vector<Vector3> &P, int n)
{
/** Create a best fit Polynomial of degree n to a set of data points.
 *
 * This function finds the best fit polynomial (in a linear least squares sense)
 * to a set of data points (x,y). The polynomial returned is of degree n.
 *
 * Inputs:
 * P     -- a vector of Vector3s: data points for approximation
 * n     -- terms in polynomial (order = n-1)
 *
 * Outputs:
 * xpoly -- a polynomial evaluated in terms of x.
 **/

    Valmatrix A(P.size(), n, 0.0);
    vector<double> b(P.size());
    vector<double> x(n);
    double xval = 0.0;

    for( size_t i = 0; i < P.size(); ++i ) {
	xval = P[i].x;
	for(int j = 0; j < n; ++j ) {
	    A.set(i, j, pow(xval, j));
	}
	b[i] = P[i].y;
    }

    least_squares_solve(A, x, b);

    vector<double> B;
   
    for( size_t i = 0; i < x.size(); ++i ) {
	B.push_back(x[i]);
    }

    return XPoly(P.front().x, P.back().x, B);
}

Best_fit_Bezier::
Best_fit_Bezier(vector<Vector3> points, int order,
		double dydxA, double dydxB, double z_pos)
    : MultivariateFunction(), points_(points),
      order_(order), dydxA_(dydxA), dydxB_(dydxB), z_pos_(z_pos) {}

Best_fit_Bezier::
~Best_fit_Bezier() {}

double
Best_fit_Bezier::
eval(vector<double> &x)
{
    const int nsamples = 1000;
    const double dt = 1.0/(nsamples-1);

    vector<Vector3> B;
    B.resize(order_+1);
    B[0] = points_[0];
    B[0].z = z_pos_;
    B.back() = points_.back();
    B.back().z = z_pos_;

    B[1].x = x[0];
    B[1].y = dydxA_*(x[0] - B[0].x) + B[0].y;
    B[1].z = z_pos_;

    B[B.size()-2].x = x.back();
    B[B.size()-2].y = dydxB_*(x.back() - B.back().x) + B.back().y;
    B[B.size()-2].z = z_pos_;

    for ( size_t i = 2; i < B.size()-2; ++i ) {
	B[i].x = x[2*(i-1)-1];
	B[i].y = x[2*(i-1)];
	B[i].z = z_pos_;
    }

    Bezier bez(B);

    vector<double> errors;
    errors.resize(points_.size()-2);
    double t = 0.0;
    Vector3 test;
    test = bez.eval(t);

    for ( size_t i = 0; i < errors.size(); ++i )
	errors[i] = vabs(points_[i+1] - test);

    double L = 0.0;
    for ( t = dt; t <= 1.0; t += dt ) {
	test = bez.eval(t);
	for ( size_t i = 0; i < errors.size(); ++i ) {
	    L = vabs(points_[i+1] - test);
	    if ( L < errors[i] )
		errors[i] = L;
	}
    }

    double sum = 0.0;
    for ( size_t i = 0; i < errors.size(); ++i )
	sum += errors[i];
    return sum;
}

Best_fit_Bezier_YZ::
Best_fit_Bezier_YZ(vector<Vector3> points, int order,
		   double dydzA, double dydzB, double x_pos)
    : MultivariateFunction(), points_(points),
      order_(order), dydzA_(dydzA), dydzB_(dydzB), x_pos_(x_pos) {}

Best_fit_Bezier_YZ::
~Best_fit_Bezier_YZ() {}

double
Best_fit_Bezier_YZ::
eval(vector<double> &z)
{
    const int nsamples = 1000;
    const double dt = 1.0/(nsamples-1);

    vector<Vector3> B;
    B.resize(order_+1);
    B[0] = points_[0];
    B[0].x = x_pos_;
    B.back() = points_.back();
    B.back().x = x_pos_;

    B[1].z = z[0];
    B[1].y = dydzA_*(z[0] - B[0].z) + B[0].y;
    B[1].x = x_pos_;

    B[B.size()-2].z = z.back();
    B[B.size()-2].y = dydzB_*(z.back() - B.back().z) + B.back().y;
    B[B.size()-2].x = x_pos_;

    for ( size_t i = 2; i < B.size()-2; ++i ) {
	B[i].z = z[2*(i-1)-1];
	B[i].y = z[2*(i-1)];
	B[i].x = x_pos_;
    }

    Bezier bez(B);

    vector<double> errors;
    errors.resize(points_.size()-2);
    double t = 0.0;
    Vector3 test;
    test = bez.eval(t);

    for ( size_t i = 0; i < errors.size(); ++i )
	errors[i] = vabs(points_[i+1] - test);

    double L = 0.0;
    for ( t = dt; t <= 1.0; t += dt ) {
	test = bez.eval(t);
	for ( size_t i = 0; i < errors.size(); ++i ) {
	    L = vabs(points_[i+1] - test);
	    if ( L < errors[i] )
		errors[i] = L;
	}
    }

    double sum = 0.0;
    for ( size_t i = 0; i < errors.size(); ++i )
	sum += errors[i];
    return sum;
}


Best_fit_Bezier3D::
Best_fit_Bezier3D(vector<Vector3> points, int order,
		  double dydxA, double dydxB,
		  double dzdxA, double dzdxB)
    : MultivariateFunction(), points_(points),
      order_(order), dydxA_(dydxA), dydxB_(dydxB),
      dzdxA_(dzdxA), dzdxB_(dzdxB) {}

Best_fit_Bezier3D::
~Best_fit_Bezier3D() {}

double
Best_fit_Bezier3D::
eval(vector<double> &x)
{
    const int nsamples = 1000;
    const double dt = 1.0/(nsamples-1);

    vector<Vector3> B;
    B.resize(order_+1);
    B[0] = points_[0];
    B.back() = points_.back();

    B[1].x = x[0];
    B[1].y = dydxA_*(x[0] - B[0].x) + B[0].y;
    B[1].z = dzdxA_*(x[0] - B[0].x) + B[0].z;

    B[B.size()-2].x = x.back();
    B[B.size()-2].y = dydxB_*(x.back() - B.back().x) + B.back().y;
    B[B.size()-2].z = dzdxB_*(x.back() - B.back().x) + B.back().z;

    for ( size_t i = 2; i < B.size()-2; ++i ) {
	B[i].x = x[3*(i-1)-2];
	B[i].y = x[3*(i-1)-1];
	B[i].z = x[3*(i-1)];
    }

    Bezier bez(B);

    vector<double> errors;
    errors.resize(points_.size()-2);
    double t = 0.0;
    Vector3 test;
    test = bez.eval(t);

    for ( size_t i = 0; i < errors.size(); ++i )
	errors[i] = vabs(points_[i+1] - test);

    double L = 0.0;
    for ( t = dt; t <= 1.0; t += dt ) {
	test = bez.eval(t);
	for ( size_t i = 0; i < errors.size(); ++i ) {
	    L = vabs(points_[i+1] - test);
	    if ( L < errors[i] )
		errors[i] = L;
	}
    }

    double sum = 0.0;
    for ( size_t i = 0; i < errors.size(); ++i )
	sum += errors[i];
    return sum;
}


int intersect_3D_lines(const Vector3 &P0, const Vector3 &T0,
		       const Vector3 &P2, const Vector3 &T2,
		       double &alf1, double &alf2,
		       Vector3 &P1)
{
    const double SMALL_VAL = 1.0e-12;
    const double TOL = 1.0e-6; // tolerance for determining if two
                               // Cartesian points are equal

    Vector3 w = P0 - P2;
    double a = dot(T0, T0);
    double b = dot(T0, T2);
    double c = dot(T2, T2);
    double d = dot(T0, w);
    double e = dot(T2, w);
    double D = a*c - b*b;
    double sc = 0.0;
    double tc = 0.0;

    int status = 0;

    if ( D < SMALL_VAL ) {
	sc = 0.0;
	tc = ( b > c ? d/b : e/c );
	return 1; // parallel or NOT coplanar
    }
    else {
	sc = (b*e - c*d) / D;
	tc = (a*e - b*d) / D;
    }
    
    alf1 = sc;
    alf2 = tc;
    P1 = P0 + sc*T0;
    Vector3 A = P2 + tc*T2;

    if ( ! equal(P1, A, TOL) ) {
	status = 1;
    }

    return status;
}
