// Author: Rowan J. Gollan
// Date: 19-May-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <valarray>
#include <fstream>

#include "../../util/source/useful.h"
#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "gpath_utils.hh"
#include "nurbs_utils.hh"
#include "surface.hh"

using namespace std;


//****f*
int curve_knot_ins(int np, int p, const std::vector<double> &UP,
		   const std::vector<Mapped_point> &Pw, double u, 
		   int k, int s, int r, int &nq, 
		   std::vector<double> &UQ, std::vector<Mapped_point> &Qw)
{
/** Compute new curve from knot insertion.
 *
 * This is Algorithm A5.1 in Piegl and Tiller (1997).  See Section 5.2.
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 * Input:
 * np    -- index of last control point in original curve
 * p     -- degree of original curve
 * UP    -- knot vector for original curve
 * Pw    -- control points (in 4D/hyperspace) for original curve
 * u     -- knot to insert
 * k     -- index for knot insertion : u into [u_k, u_k+1)
 * s     -- multiplicity of u in original UP
 * r     -- no. times to insert u
 *
 * Output:
 * nq    -- no. control points in new curve
 * UQ    -- knot vector for new curve
 * Qw    -- control points (in 4D/hyperspace) for new curve
 *
 * Return:
 * SUCCESS <or>
 * FAILURE
 *
 **/
    int mp  = np + p + 1;
    nq = np + r;
    Qw.resize(Pw.size()+r);
    UQ.resize(UP.size()+r);
    vector<Mapped_point> Rw(p+1);
    
    // Load new knot vector
    for ( int i = 0; i <= k; ++i ) UQ[i] = UP[i];
    for ( int i = 1; i <= r; ++i ) UQ[k+1] = u;
    for ( int i = k+1; i <= mp; ++i ) UQ[i+r] = UP[i];

    // Save unaltered control points
    for ( int i = 0; i <= k; ++i ) Qw[i] = Pw[i];
    for ( int i = k-s; i <= np; ++i ) Qw[i+r] = Pw[i];
    for ( int i = 0; i <= p-s; ++i ) Rw[i] = Pw[k-p+i];

    // Insert the knot 'r' times
    int L = 0;
    for ( int j = 1; j <= r; ++j ) {
	L = k - p + j;
	for ( int i = 0; i <= p-j-s; ++i ) {
	    double alpha = (u - UP[L+i])/(UP[i+k+1]-UP[L+i]);
	    Rw[i] = alpha*Rw[i+1] + (1.0 - alpha)*Rw[i];
	}
	Qw[L] = Rw[0];
	Qw[k+r-j-s] = Rw[p-j-s];
    }

    for ( int i = L+1; i < k-s; ++i ) {
	Qw[i] = Rw[i-L];
    }

    return SUCCESS;

}

//****f*
int refine_knot_vector_curve(int n, int p, const vector<double> &U,
			     const vector<Mapped_point> &Pw,
			     const vector<double> &X, int r,
			     vector<double> &Ubar,
			     vector<Mapped_point> &Qw)
{
/** Refine a curve knot vector. 
 *
 * This appears as Algorithm A5.4 in Piegl and Tiller (1997).
 * They reference Boehm and Prautzsch (1985) for this algorithm.
 *
 * References:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 * Boehm, W. and Prautzsch, H. (1985)
 * The insertion algorithm.
 * CAD vol. 17, no. 2, pp. 58--59
 *
 * Input:
 * n     -- index of last control point in original curve
 * p     -- degree of curve
 * U     -- original knot vector
 * Pw    -- control points (in 4D/hyperspace) for original curve
 * X     -- list of knots to be inserted
 * r     -- index of last knot to be inserted
 *
 * Output:
 * Ubar  -- new knot vector, after insertion
 * Qw    -- control points (in 4D/hyperspace) for new curve
 *
 **/

//    cout << "--- knot_refinement called ---\n";
    const double ZERO_TOL = 1.0e-15;

    int m = n + p + 1;
    int a = find_span(X[0], n, p, U);
    int b = find_span(X[r], n, p, U);
    b = b + 1;
    
    Ubar.resize(U.size()+X.size());
    Qw.resize(Pw.size()+X.size());

    for ( int j = 0; j <= a-p; ++j ) Qw[j] = Pw[j];
    for ( int j = b-1; j <= n; ++j ) Qw[j+r+1] = Pw[j];
    for ( int j = 0; j <= a; ++j ) Ubar[j] = U[j];
    for ( int j = b+p; j <= m; ++j ) Ubar[j+r+1] = U[j];

    int i = b+p-1;
    int k = b+p+r;

//     cout << "Pw:\n";
//     for ( size_t g = 0; g < Pw.size(); ++g ) {
// 	cout << "i= " << g << " w= " << Pw[g].w << endl;
//     }

//     cout << "Qw:\n";
//     for ( size_t g = 0; g < Qw.size(); ++g ) {
// 	cout << "i= " << g << " w= " << Qw[g].w << endl;
//     }

    
    for ( int j = r; j >= 0; --j ) {
	
	while ( X[j] <= U[i] && i > a ) {
	    Qw[k-p-1] = Pw[i-p-1];
	    Ubar[k] = U[i];
	    k = k - 1;
	    i = i - 1;
	}

	Qw[k-p-1] = Qw[k-p];

	for ( int l = 1; l <= p; ++l ) {
	    int ind = k-p+l;
	    double alfa = Ubar[k+l] - X[j];
	    if ( fabs(alfa) < ZERO_TOL ) {
		Qw[ind-1] = Qw[ind];
	    }
	    else {
		alfa = alfa / (Ubar[k+l] - U[i-p+l]);
		Qw[ind-1] = alfa*Qw[ind-1] + (1.0 - alfa)*Qw[ind];
	    }
	}
	Ubar[k] = X[j];
	k = k - 1;
    }

    return SUCCESS;
}

//****f*
int estimate_tangents(const vector<Vector3> &Q, vector<Vector3> &T,
		      bool preserve_corners)
{
/** Estimate the tangents for a set of points.
 *
 * The tangent estimate is based on Eqn 9.29 and 9.31 in
 * Piegl and Tiller (1997).  See section 9.3.1.
 *
 * Input:
 * Q      -- vector of points to interpolate
 * T      -- vector to hold tangents
 * preserve_corners -- flag to indicate if corners should be preserved
 *
 * Output:
 * T      -- is filled with tangents
 *
 * Return:
 * SUCCESS -- if successful, else
 * FAILURE
 *
 **/

    const double ZERO_TOL = 1.0e-12;

    vector<Vector3> q(Q.size()+3);

    for ( size_t k = 1; k < Q.size(); ++k ) {
	// The first index is offset by one from the algorithm
	q[k+1] = Q[k] - Q[k-1];
    }
    // Special cases
    q[1] = 2*q[2] - q[3];
    q[0] = 2*q[1] - q[2];
    q[q.size()-2] = 2*q[q.size()-3] - q[q.size()-4];
    q[q.size()-1] = 2*q[q.size()-2] - q[q.size()-3];
    
    double alpha_k;
    double numer;
    double denom;
    for ( size_t k = 0; k < T.size(); ++k ) {
	// All indices for q[k] are increased by one because of the offset
	// compared to the equations in Piegl and Tiller.
	// If C++ let you start an array from -1, then we'd be set, but alas...
	numer = vabs(cross(q[k], q[k+1]));
	denom = vabs(cross(q[k], q[k+1])) + vabs(cross(q[k+2], q[k+3]));

	if ( denom <= ZERO_TOL ) {
	    if ( preserve_corners ) {
		alpha_k = 1.0;
	    }
	    else {
		alpha_k = 0.5;
	    }
	}
	else {
	    alpha_k = numer/denom;
	}

	T[k] = unit((1.0 - alpha_k)*q[k+1] + alpha_k*q[k+2]);
    }

    return SUCCESS;

}

//****f*
Nurbs local_rat_quad_curve_interp(const vector<Vector3> &Q, bool preserve_corners)
{
/** Build a Nurbs curve that interpolates data points.
 *
 * Each segment is a rational quadratic curve (or curves, if necessary)
 * which pass through the end points. The segments maintain G1 continuity,
 * that is, they appear visually smooth.
 *
 * See Sections 9.3.2 and 9.3.3 of Piegl and Tiller (1997).
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 * Input:
 * Q      -- vector of points to interpolate
 * preserve_corners -- flag to indicate if corners should be preserved.
 *
 * Output:
 * n      -- a Nurbs curve interpolating the points
 *
 **/

    vector<Vector3> T(Q.size());
    if ( estimate_tangents(Q, T, preserve_corners) != SUCCESS ) {
	cout << "local_rat_quad_curve_interp():\n";
	cout << "There was a problem estimating the tangent vectors\n";
	cout << "for the set of given points.\n";
	cout << "Bailing out!\n";
	exit(NUMERICAL_ERROR);
    }

    vector<Vector3> C;
    C.push_back(Q[0]);
    
    vector<double> weights;
    weights.push_back(1.0);

    vector<double> U;
    U.push_back(0.0); U.push_back(0.0); U.push_back(0.0);

    vector<Vector3> R;
    vector<double> w;

    double chord = 0.0;

    for ( size_t k = 1; k < Q.size(); ++k ) {

	int result = quad_curve_interp(Q[k-1], T[k-1], Q[k], T[k],
				       R, w);

	if ( result != 1 ) {
	    cout << "local_rat_quad_curve_interp():\n";
	    cout << "NUMBER OF SEGMENTS NOT EQUAL TO 1.\n";
	    cout << "NOT HANDLED CORRECTLY YET.\n";
	    cout << "Bailing out!\n";
	    exit(1);
	}

	C.push_back(R[0]);
	weights.push_back(w[0]);

	C.push_back(Q[k]);
	weights.push_back(1.0);

	chord += vabs(Q[k] - Q[k-1]);

	U.push_back(chord); U.push_back(chord);
    }

    // Need one more knot at end.
    U.push_back(chord);

    // Normalise knot vector.
    for ( size_t i = 3; i < U.size(); ++i ) {
	U[i] = U[i]/U.back();
    }

    // Create Nurbs curve
    
    int p = 2;
    return Nurbs(C, weights, p, U);
}

//****f*
int quad_curve_interp(const Vector3 &Q0, const Vector3 &T0,
		      const Vector3 &Q1, const Vector3 &T1,
		      vector<Vector3> &R,
		      vector<double> &w)
{
/** Find a local rational quadratic curve interpolation to Q0 and Q1.
 *
 * This function returns middle control point and weight for a rational
 * quadratic curve interpolating Q0 and Q1.  When the tangents are parallel,
 * but not parallel to the chord Q0-Q1, the interpolating curve is split
 * into 2 segments, each rational quadratic curves.
 *
 * This implementation follows the description that appears in 
 * Sections 9.3.2 and 9.3.3 of Piegl and Tiller.
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 * Input:
 * Q0    -- start point
 * T0    -- unit tangent vector at start point
 * Q1    -- end point
 * T1    -- unit tangent vector at end point
 *
 * Output:
 * R     -- mid control point (or points when extra segments present)
 * w     -- weight at mid control point (or weights...)
 *
 * Returns:
 * -1    -- for an error
 * 1     -- when single segment computed
 * 2     -- when two segments computed
 *
 **/

    const double TOL = 1.0e-12;

    if ( fabs(vabs(T0) - 1.0) > 1.0e-12 ) {
	cout << "quad_curve_interp():\n";
	cout << "T0 in NOT a unit vector: " << T0 << endl;
	cout << "|T0| = " << vabs(T0) << endl;
	cout << "A unit tangent vector is required.\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    if ( fabs(vabs(T1) - 1.0) > 1.0e-12 ) {
	cout << "quad_curve_interp():\n";
	cout << "T1 in NOT a unit vector: " << T1 << endl;
	cout << "|T1| = " << vabs(T1) << endl;
	cout << "A unit tangent vector is required.\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    // To be sure, clear the passed in vectors
    R.clear();
    w.clear();

    double alf1 = 0.0;
    double alf2 = 0.0;
    Vector3 Rk;
    double wk;

    int intersect = intersect_3D_lines(Q0, T0, Q1, T1, alf1, alf2, Rk);

    // Special case for corners.
    if ( intersect == 0 && 
	 ( equal(Q0, Rk, TOL) || equal(Q1, Rk, TOL) ) ) {
	// This is the same as Case 2.
	Rk = 0.5*(Q0 + Q1);

	if ( quad_curve_interp_weight(Q0, Rk, Q1, wk) != SUCCESS ) {
// 	    cout << "quad_curve_interp():\n";
// 	    cout << "Problem computing weight at a corner.\n";
	    return -1;
	}

	R.push_back(Rk);
	w.push_back(wk);

	return 1;
    }
    

    // Case 1. Intersection exists, and
    //         alf1 > 0, alf2 < 0

    if ( intersect == 0  &&
	 alf1 > 0.0 &&
	 alf2 < 0.0 ) {

	if ( quad_curve_interp_weight(Q0, Rk, Q1, wk) != SUCCESS ) {
// 	    cout << "quad_curve_interp():\n";
// 	    cout << "Problem computing weight when intersection exists.\n";
	    return -1;
	}

	R.push_back(Rk);
	w.push_back(wk);

	return 1;
    }

    // Case 2. T0 and T1 are parallel
    Vector3 Q0Q1 = Q1 - Q0;

    if ( equal(unit(Q0Q1), T0) &&
	 equal(unit(Q0Q1), T1) ) {
	Rk = 0.5*(Q0 + Q1);

	if ( quad_curve_interp_weight(Q0, Rk, Q1, wk) != SUCCESS ) {
// 	    cout << "quad_curve_interp():\n";
// 	    cout << "Problem computing weight when end slopes are parallel.\n";
	    return -1;
	}

	R.push_back(Rk);
	w.push_back(wk);

	return 1;
    }

    double gamma1, gamma2;
    if ( intersect != 0 ) { // lines don't intersect.
	gamma1 = 0.5*vabs(Q0Q1);
	gamma2 = gamma1;
    }
    else {
	Vector3 Q0Rk = Rk - Q0;
	Vector3 Q1Rk = Rk - Q1;

	double cos_theta0 = dot(Q0Rk, Q0Q1)/(vabs(Q0Rk)*vabs(Q0Q1));
	double cos_theta1 = dot(Q1Rk, Q0Q1)/(vabs(Q1Rk)*vabs(Q0Q1));

	double alpha = 2.0/3.0; // As per p. 392 of Piegl and Tiller
	gamma1 = 0.5*(vabs(Q0Q1))/(1.0 + alpha*cos_theta1 + (1.0 - alpha)*cos_theta0);
	gamma2 = 0.5*(vabs(Q0Q1))/(1.0 + alpha*cos_theta0 + (1.0 - alpha)*cos_theta1);
    }

    Vector3 R1 = Q0 + gamma1*T0;
    Vector3 R2 = Q1 - gamma2*T1;
    Vector3 Qk = (gamma1*R2 + gamma2*R1)/(gamma1 + gamma2);

    double w1;
    if ( quad_curve_interp_weight(Q0, R1, Qk, w1) != SUCCESS ) {
// 	cout << "quad_curve_interp_weight():\n";
// 	cout << "Problem computing weight 1.\n";
	return -1;
    }

    double w2;
    if ( quad_curve_interp_weight(Qk, R2, Q1, w2) != SUCCESS ) {
// 	cout << "quad_curve_interp_weight():\n";
// 	cout << "Problem computing weight 2.\n";
	return -1;
    }

    R.push_back(R1); R.push_back(Qk); R.push_back(R2);
    w.push_back(w1); w.push_back(1.0); w.push_back(w2);

    return 2;  // 2 segments computed

}

int quad_curve_interp_weight(const Vector3 &Q0, const Vector3 &R1, const Vector3 &Q1,
			     double &w)
{
    const double TOL = 1.0e-9; // for checking isosceles triangle

    Vector3 Q0R1 = R1 - Q0;
    Vector3 Q1R1 = R1 - Q1;
    Vector3 Q0Q1 = Q1 - Q0;

    // Case 1. Q0, R1 and Q1 are collinear
    if ( equal(unit(Q0R1), unit(Q1R1)) ) {
	w = 1.0;
	return SUCCESS;
    }

    // Case 2.
    if ( fabs(vabs(Q0R1) - vabs(Q1R1)) < TOL ) {
	w = dot(Q0R1, Q0Q1) / (vabs(Q0R1) * vabs(Q0Q1));
	return SUCCESS;
    }

    // Case 3.
    // a. 
    Vector3 M = 0.5*(Q0 + Q1);
    
    // b.
    Vector3 MR1 = R1 - M;
    double ratio = vabs(Q0R1)/vabs(Q0Q1);
    double frac = 1.0/(ratio + 1.0);
    Vector3 D = Q1 + frac*Q1R1;
    Vector3 Q0D = D - Q0;
    Vector3 S1;
    double alf1 = 0.0;
    double alf2 = 0.0;

    if ( intersect_3D_lines(Q0, unit(Q0D), M, unit(MR1), alf1, alf2, S1) != 0 ) {
// 	cout << "quad_curve_interp_weight():\n";
// 	cout << "Problem finding intersection with bisector 1.\n";
// 	cout << "Q0= " << Q0 << endl;
// 	cout << "T0= " << Q0D << endl;
// 	cout << "Q1= " << M << endl;
// 	cout << "T1= " << MR1 << endl;
	return FAILURE;
    }
    
    // c.
    ratio = vabs(Q1R1)/vabs(Q0Q1);
    frac = 1.0/(ratio + 1.0);
    D = Q0 + frac*Q0R1;
    Vector3 Q1D = D - Q1;
    Vector3 S2;
    if ( intersect_3D_lines(Q1, unit(Q1D), M, unit(MR1), alf1, alf2, S2) != 0 ) {
// 	cout << "quad_curve_interp_weight():\n";
// 	cout << "Problem finding intersection with bisector 2.\n";
// 	cout << "Q0= " << Q1 << endl;
// 	cout << "T0= " << unit(Q1D) << endl;
// 	cout << "Q1= " << M << endl;
// 	cout << "T1= " << unit(MR1) << endl;
	return FAILURE;
    }
    
    // d.
    Vector3 S = 0.5*(S1 + S2);
    Vector3 MS = S - M;
    double s = vabs(MS)/vabs(MR1);
    w = s / (1.0 - s);

    return SUCCESS;
}

//****f*
int make_one_arc(const Vector3 &P0, const Vector3 &T0, const Vector3 &P2, const Vector3 &T2,
		 const Vector3 &P, Vector3 &P1, double &w1)
{
/** Make an arc leaving P0 at slope T0 and joining P2 at slope T2, passing through P.
 *
 * This function returns middle control point and weight for a rational
 * quadratic Bezier representing an arc.
 *
 * Input:
 * P0    -- start point
 * T0    -- unit tangent vector at start point
 * P2    -- end point
 * T2    -- unit tangent vector at end point
 * P     -- interior point the arc passes through
 *
 * Output:
 * P1    -- mid control point
 * w1    -- weight at mid control point
 *
 * Returns:
 * SUCCESS or
 * FAIL
 * 
 * This is an implementation of Algorithm A7.2 (p. 314) from Piegl and Tiller.
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 **/

    Vector3 V02 = P2 - P0;
    double dummy = 0.0;
    Vector3 dummy_P;
    double alf0 = 0.0;
    double alf2 = 0.0;
    int result = intersect_3D_lines(P0, T0, P2, T2, dummy, dummy, P1);

    if ( result == 0 ) {
	// finite control point
	Vector3 V1P = P - P1;
	intersect_3D_lines(P1, V1P, P0, V02, alf0, alf2, dummy_P);
	double a = sqrt(alf2/(1.0 - alf2));
	double u = a/(1.0 + a);
	double numer = (1.0 - u)*(1.0 - u)*dot(P-P0, P1-P) + u*u*dot(P-P2, P1-P);
	double denom = 2.0*u*(1.0 - u)*dot(P1-P, P1-P);
	w1 = numer/denom;
	if ( std::isnan(w1) || std::isinf(w1) ) {
	    return FAILURE;
	}
    }
    else {
	// infinite control point
	w1 = 0.0;
	intersect_3D_lines(P, T0, P0, V02, alf0, alf2, dummy_P);
	double a = sqrt(alf2/(1.0 - alf2));
	double u = a/(1.0 + a);
	double b = 2.0*u*(1.0-u);
	b = -alf0*(1.0 - b)/b;
	P1 = b*T0;
    }

    return SUCCESS;

}

//****f*
void deriv_basis_funs(int i, double u, int p, int n,
		      const vector<double> &U, vector<vector<double> > &ders)
{
/** Compute nonzero basis function and their derivatives.
 *
 * Algorithm A2.3
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 **/

    vector<vector<double> > ndu;
    ndu.resize(p+1);
    for ( int j = 0; j < p+1; ++j ) ndu[j].resize(p+1);
    ndu[0][0] = 1.0;

    vector<vector<double> > a;
    a.resize(2);
    a[0].resize(p+1); a[1].resize(p+1);

    vector<double> left(p+1, 0.0);
    vector<double> right(p+1, 0.0);

    for ( int j = 1; j <= p; ++j ) {
	left[j] = u - U[i+1-j];
	right[j] = U[i+j] - u;

	double saved = 0.0;

	for ( int r = 0; r < j; ++r ) {
	    ndu[j][r] = right[r+1] + left[j-r];
	    double temp = ndu[r][j-1]/ndu[j][r];
	    
	    ndu[r][j] = saved + right[r+1]*temp;
	    saved = left[j-r]*temp;
	}

	ndu[j][j] = saved;
    }

    for ( int j = 0; j <= p; ++j )
	ders[0][j] = ndu[j][p];

    for ( int r = 0; r <= p; ++r ) {
	int s1 = 0;
	int s2 = 1;
	a[0][0] = 1.0;

	for ( int k = 1; k <= n; ++k ) {
	    double d = 0.0;
	    int rk = r - k;
	    int pk = p - k;
	    
	    if ( r >= k ) {
		a[s2][0] = a[s1][0]/ndu[pk+1][rk];
		d = a[s2][0]*ndu[rk][pk];
	    }
	    
	    int j1, j2;
	    if ( rk >= -1 )
		j1 = 1;
	    else
		j1 = -rk;
	    
	    if ( r-1 <= pk )
		j2 = k - 1;
	    else
		j2 = p - r;

	    int j;
	    for ( j = j1; j <= j2; ++j ) {
		a[s2][j] = (a[s1][j] - a[s1][j-1])/(ndu[pk+1][rk+j]);
		d += a[s2][j]*ndu[rk+j][pk];
	    }

	    if ( r <= pk ) {
		a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
		d += a[s2][k]*ndu[rk+j][pk];
	    }

	    ders[k][r] = d;

	    j = s1; s1 = s2; s2 = j;
	}
    }

    int r = p;

    for ( int k = 1; k <= n; ++k ) {
	for ( int j = 0; j <= p; ++j )
	    ders[k][j] *= r;

	r *= (p - k);
    }
}

//****f*
void curve_derivs(int n, int p, const vector<double> &U,
		  const vector<Mapped_point> &Pw,
		  double u, int d,
		  vector<Mapped_point> &CKw)
{
/** Compute nonzero basis function and their derivatives.
 *
 * Algorithm A3.2
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 **/

    int du = min(d, p);
    for ( int k = p+1; k <= d; ++k ) {
	CKw[k] = CKw[k] * 0.0;
    }

    int span = find_span(u, n, p, U);
    
    vector<vector<double> > nders;
    nders.resize(p+1);
    for ( size_t j = 0; j < nders.size(); ++j ) nders[j].resize(p+1);

    deriv_basis_funs(span, u, p, du, U, nders);
    
    for ( int k = 0; k <= du; ++k ) {

	CKw[k] = CKw[k] * 0.0;

	for ( int j = 0; j <= p; ++j ) {
	    CKw[k] = CKw[k] + nders[k][j]*Pw[span-p+j];
	}
    }
}

int factorial(int n)
{
    if ( n <= 1 )
	return 1;
    else
	return n*factorial(n - 1);
}

int bin_coeff(int n, int k)
{
    return factorial(n) / (factorial(k) * factorial(n - k));
}

//****f*
void rat_curve_derivs(const vector<Vector3> &Aders, const vector<double> &wders,
		      int d, vector<Vector3> &CK)
{
/** Compute C(u) derivatives based on Cw(u) derivatives.
 *
 * Algorithm A4.2
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 **/
    for ( int k = 0; k <= d; ++k ) {

	Vector3 v = Aders[k];
	
	for ( int i = 1; i <= k; ++i ) 
	    v = v - bin_coeff(k, i)*wders[i]*CK[k-i];
	
	CK[k] = v/wders[0];
    }
}

void rat_curve_derivs(const vector<Mapped_point> &CKw, int d, vector<Vector3> &CK)
{
    vector<Vector3> Aders;
    Aders.resize(CKw.size());
    
    vector<double> wders;
    wders.resize(CKw.size());

    for ( size_t i = 0; i < CKw.size(); ++i ) {
	Aders[i].x = CKw[i].wx; Aders[i].y = CKw[i].wy; Aders[i].z = CKw[i].wz;
	wders[i] = CKw[i].w;
    }

    rat_curve_derivs(Aders, wders, d, CK);

}

double dist_point_projection(const Vector3 &P, const Nurbs &C,
			     double &t_found, Vector3 &C_found,
			     double eps1, double eps2)
{

    const int SAMPLE_NO = 20;
    const int MAX_ITERATIONS = 100;

    // 0. Sample some points to find a good starting u
    double t_min = 0.0;
    Vector3 P_min = C.eval(t_min);

    for ( int i = 1; i < SAMPLE_NO; ++i ) {
	double t_test = i / (SAMPLE_NO - 1.0);
	Vector3 P_test = C.eval(t_test);

	if ( vabs(P_test - P) < vabs(P_min - P) ) {
	    t_min = t_test;
	    P_min = P_test;
	}
    }

    // 1. Begin iterating starting with t_min
    double t0 = t_min;
    double t1;
    int iter;

    bool criterion1 = false;
    bool criterion2 = false;
    bool criterion4 = false;
	
    Vector3 C0; // value at t0
    Vector3 C1; // first deriv
    Vector3 C2; // second deriv
    Vector3 PC0; // C0 - P
    double dot_C1_PC0;
    double dot_C2_PC0;
    double abs_PC0;
    double abs_C1;
    double cos_T;
    double u0, u1;

    for ( iter = 0; iter < MAX_ITERATIONS; ++iter ) {

	C0 = C.eval(t0);
	C1 = C.dpdu(t0);
	C2 = C.d2pdu2(t0);
	PC0 = C0 - P;
	dot_C1_PC0 = dot(C1, PC0);
	dot_C2_PC0 = dot(C2, PC0);
	abs_PC0 = vabs(PC0);
	abs_C1 = vabs(C1);

	t1 = t0 - dot_C1_PC0 / ( dot_C2_PC0 + abs_C1*abs_C1);
	
	// Test criterion 1.
	if ( abs_PC0 <= eps1 )
	    criterion1 = true;

	// Test criterion 2
	cos_T = fabs(dot_C1_PC0) / (abs_C1 * abs_PC0);
	if ( cos_T <= eps2 )
	    criterion2 = true;

	if ( criterion1 && criterion2 ) {
	    t_found = t0; C_found = C0;
	    return abs_PC0;
	}
	// else go on..

	// Test criterion 3.
	// ASSUMPTION: Our curves are NOT closed.
	// Need to fix this if that assumption changes.
	//
	// ASSUMPTION: Our curve extends from t-0.0 --> t=1.0
	//

	if ( t1 < 0.0 )
	    t1 = 0.0;

	if ( t1 > 1.0 )
	    t1 = 1.0;

	// Test criterion 4.
	u0 = C.map_t2u(t0);
	u1 = C.map_t2u(t1);

	if ( vabs((u1 - u0)*C1) <= eps1 )
	    criterion4 = true;

	if ( criterion1 || criterion2 || criterion4 ) {
	    t_found = t0; C_found = C0;
	    return abs_PC0;
	}

	t0 = t1;
	criterion1 = false;
	criterion2 = false;
	criterion4 = false;
    }

    // If we get out here, we did not converge.
    // Eventually, handle gracefully.
    // For the moment.. just crash.

    cout << "dist_point_projection():\n";
    cout << "Iterations did not converge.\n";
    cout << "Bailing out!\n";
    exit(1);
}

//****f*
int fit_with_conic(int ks, int ke, const vector<Vector3> &Q,
		   const Vector3 &Ts, const Vector3 &Te,
		   double E, Mapped_point &Rw)
{
/** Fit to tolerance E with conic segment.
 *
 * Algorithm A9.11
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 **/

    // Not really sure what are reasonable values for these
    // Piegl and Tiller give no clues.
    const double WMIN = 0.001;
    const double WMAX = 10.0; 

    if ( (ke - ks) == 1 ) {
	// No interior segments.
	// Fit an interpolating segment.
	vector<Vector3> R;
	vector<double> w;

	int result = quad_curve_interp(Q[ks], Ts, Q[ke], Te,
				       R, w);

	if ( result == 1 ) {
	    Rw.wx = w[0]*R[0].x;
	    Rw.wy = w[0]*R[0].y;
	    Rw.wz = w[0]*R[0].z;
	    Rw.w = w[0];
	    return 1;
	}
	else {
	    // Try just using a straignt line segment.
	    // Do this by setting the tangents along the line
	    // joining Q[ks] and Q[ke].
	    Vector3 T = unit(Q[ke] - Q[ks]);
	    result = quad_curve_interp(Q[ks], T, Q[ke], T,
				       R, w);
	    if ( result != 1 ) {
		cout << "fit_with_conic():\n";
		cout << "No simple interpolation found for two points\n";
		cout << "without any interior points.\n";
		cout << "Failed at fitting a straight line also.\n";
		cout << "Q0= " << Q[ks] << " T0= " << Ts << endl;
		cout << "Q1= " << Q[ke] << " Te= " << Te << endl;
		cout << "Bailing out!\n";
		exit(1);
	    }
	    Rw.wx = w[0]*R[0].x;
	    Rw.wy = w[0]*R[0].y;
	    Rw.wz = w[0]*R[0].z;
	    Rw.w = w[0];
	    return 1;
	}
    }

    double alf1, alf2;
    Vector3 R;
    int i = intersect_3D_lines(Q[ks], Ts, Q[ke], Te, alf1, alf2, R);

    if ( i != 0 ) {
	// No intersection
	// Test for collinearity of points
	Vector3 Q0Q1 = Q[ke] - Q[ks];

	for ( int k = ks+1; k <= ke-1; ++k ) {
	    if ( ! equal(unit(Q0Q1), unit(Q[k] - Q[ks])) ) {
		// not collinear
		return 0;
	    }
	}

	R = 0.5*(Q[ks] + Q[ke]);
	Rw.wx = R.x; Rw.wy = R.y; Rw.wz = R.z; Rw.w = 1.0;
	return 1;
    }

    if ( alf1 <= 0.0 || alf2 >= 0.0 ) {
	return 0;
    }

    double s = 0.0;
    Vector3 V = Q[ke] - Q[ks];
    Vector3 dummy_P;

    for ( int k = ks+1; k <= ke-1; ++k ) {
	// Get conic which passes through each internal
	// point in turn, and form some kind of "average" weight
	Vector3 V1 = Q[k] - R;
	int j = intersect_3D_lines(Q[ks], V, R, V1, alf1, alf2, dummy_P);
	if ( j != 0 || alf1 <= 0.0 || alf1 >= 1.0 || alf2 <= 0.0 ) {
	    return 0;
	}
	
	double wk;
	make_one_arc(Q[ks], Ts, Q[ke], Te, Q[k], dummy_P, wk);

	s = s + wk/(1.0 + wk);
    }

    s = s/(ke-ks-1);
    double w = s/(1.0 - s);

    if ( w < WMIN || w > WMAX ) {
	return 0;
    }

    vector<Vector3> B(3);
    B[0] = Q[ks]; B[1]= R; B[2] = Q[ke];

    vector<double> weights(3, 1.0);
    weights[1] = w;

    vector<double> U(6, 0.0);
    U[3] = 1.0; U[4] = 1.0; U[5] = 1.0;

    int p = 2;

    Nurbs n(B, weights, p, U);

    // Check our curve is close enough
    double dist;
    double t;
    for ( int k = ks+1; k <= ke-1; ++k ) {
	dist = dist_point_projection(Q[k], n, t, dummy_P);
	if ( dist > E ) {
	    return 0;
	}
    }

    Rw.wx = w*R.x; Rw.wy = w*R.y; Rw.wz = w*R.z; Rw.w = w;
    return 1;

}

int search_max_span(int ks, int ke,
		    const vector<Vector3> &Q, const vector<Vector3> &T,
		    double E)
{

    int high = ke;
    int low = ks;
    int mid;
    Mapped_point dummy;

    if ( fit_with_conic(ks, high, Q, T[ks], T[high], E, dummy) == 1 ) {
	return high;
    }

    mid = low + (high - low)/2;

    while ( (high - low) > 1 ) {
	
	if ( fit_with_conic(ks, mid, Q, T[ks], T[mid], E, dummy) == 1 ) {
	    low = mid;
	}
	else {
	    high = mid;
	}

	mid = low + (high - low)/2;
    }

    return mid;


}

//****f*
Nurbs quad_curve_approx(const vector<Vector3> &Q, double E, int Kmax,
			bool preserve_corners)
{
/** Fit a (segmented) quadratic curve to a set of data.
 *
 * This fit approximates the data, that is, it's not guaranteed
 * to go through all of the points.  If you need to pass through
 * the points, then try an interpolation routine.
 *
 * This is the algorithm as described on p. 438 of Piegl and 
 * Tiller (1997).
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 **/

    if ( Kmax <= 0 )
	Kmax = Q.size()-1;

    const double ZERO_TOL = 1.0e-12;

    // 0. Estimate the tangent vectors at the points.
    vector<Vector3> T(Q.size());

    if ( estimate_tangents(Q, T, preserve_corners) != SUCCESS ) {
	cout << "quad_curve_approx():\n";
	cout << "There was a problem estimating the tangent vectors\n";
	cout << "for the set of given points.\n";
	cout << "Bailing out!\n";
	exit(NUMERICAL_ERROR);
    }

    // 1. Start building segments
    vector<Vector3> C; // to store control points
    C.push_back(Q[0]); // Add first point, as it's definitely there

    vector<double> weights;
    weights.push_back(1.0);

    int ks = 0;
    int m = Q.size()-1;
    int ke = min(ks + Kmax, m);
    int k_found;
    int result;
    double chord = 0.0;
    Mapped_point Rw;
    Vector3 R;
    
    vector<double> U;
    U.push_back(0.0); U.push_back(0.0); U.push_back(0.0);

    while ( ks != m ) {
	k_found = search_max_span(ks, ke, Q, T, E);
	
	result = fit_with_conic(ks, k_found, Q, T[ks], T[k_found],
				E, Rw);

	if ( result != 1 ) {
	    cout << "quad_curve_approx():\n";
	    cout << "Problem approximating span.\n";
	    cout << "ks= " << ks << " ke= " << ke << endl;
	    cout << "Bailing out!\n";
	    exit(1);
	}

	if ( Rw.w > ZERO_TOL ) {
	    R.x = Rw.wx/Rw.w; R.y = Rw.wy/Rw.w; R.z = Rw.wz/Rw.w;
	}
	else {
	    // Infinite control point.
	    R.x = Rw.wx; R.y = Rw.wy; R.z = Rw.wz;
	}

	C.push_back(R);
	weights.push_back(Rw.w);

	C.push_back(Q[k_found]);
	weights.push_back(1.0);

	chord += vabs(Q[k_found] - Q[ks]);
	U.push_back(chord); U.push_back(chord);

	ks = k_found;
	ke = min(ks + Kmax, m);
	
    }
    
    U.push_back(chord); // Need three repeated knots at end.

    // Normalise knot vector.
    for ( size_t i = 0; i < U.size(); ++i ) {
	U[i] = U[i]/U.back();
    }

    int p = 2;
    return Nurbs(C, weights, p, U);
    
}

//****f*
int interp_homogeneous_points(const vector<Mapped_point> &Qw,
			      int p, const vector<double> &u,
			      const vector<double> &U,
			      vector<Mapped_point> &Pw)
{
/** Interpolate a set of homogeneous points with a p degree curve.
 *
 * This function implements the description of an algorithm
 * given by Piegl and Tiller (1997) in Section 9.2.1 of their
 * text.  The degree of interpolating curve is chosen by
 * the caller. Also, there are several ways to get the parameters {u_k} and
 * the knot vector.  The caller needs to choose these also.
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 * Input:
 * Qw    -- vector of points in homogenous space
 * p     -- degree for interpolating curve
 * u     -- vector of parameters (corresponding to supplied points)
 * U     -- corresponding knot vector
 *
 * Output:
 * Pw    -- vector of control points for interpolating curve
 *          in homogenous space
 *
 **/

  //   for ( size_t i = 0; i < Qw.size(); ++i ) {
// 	cout << "Qw[" << i << "]: \n";
// 	cout << "wx= " << Qw[i].wx << " wy= " << Qw[i].wy << " wz= " << Qw[i].wz << " w= " << Qw[i].w << endl;
//     }
       

    Pw.resize(Qw.size());

    if ( Qw.size() != u.size() ) {
	cout << "interp_homogeneous_points():\n";
	cout << "Error: the supplied points and parameters at\n";
	cout << "which to evaluate the interpolants do not match\n";
	cout << "in size: Qw.size()= " << Qw.size() << " u.size()= " << u.size() << endl;
	cout << "Bailing out!\n";
    }

    int n = Qw.size()-1;
    int m = p + n + 1;

    if ( int(U.size()) != m+1 ) {
	cout << "interp_homogeneous_points():\n";
	cout << "Error: the supplied knot vector is not of the correct size.\n";
	cout << "Expected size: " << m+1 << " Actual size= " << U.size() << endl;
	cout << "Bailing out!\n";
    }

    // 2. Coefficient matrix
    Valmatrix A(Qw.size(), Qw.size(), 0.0);
    for ( size_t i = 0; i < Qw.size(); ++i ) {
	for ( size_t k = 0; k < Qw.size(); ++k ) {
	    A.set(i, k, one_basis_fun(p, U, k, u[i]));
	    //	    cout << "A[" << i << "][" << k << "]= " << A.get(i, k) << " ";
	}
	//	cout << endl;

    }

    // Store a copy of A because A is mangled
    // during call to solver
    Valmatrix A_store = A;

    // 3. Solve the system 4 times:
    //    1. for wx components
    //    2. for wy components
    //    3. for wz components
    //    4. for w components
    
    valarray<double> x(Qw.size()), b(Qw.size());

    // wx
    for ( size_t i = 0; i < b.size(); ++i ) b[i] = Qw[i].wx;
    gaussian_elimination(A, x, b);
    for ( size_t i = 0; i < x.size(); ++i ) {
	Pw[i].wx = x[i];
    }

    // wy
    A = A_store;
    for ( size_t i = 0; i < b.size(); ++i ) b[i] = Qw[i].wy;
    gaussian_elimination(A, x, b);
    for ( size_t i = 0; i < x.size(); ++i ) {
	Pw[i].wy = x[i];
    }

    // wz
    A = A_store;
    for ( size_t i = 0; i < b.size(); ++i ) b[i] = Qw[i].wz;
    gaussian_elimination(A, x, b);
    for ( size_t i = 0; i < x.size(); ++i ) {
	Pw[i].wz = x[i];
    }

    // w
    A = A_store;
    for ( size_t i = 0; i < b.size(); ++i ) b[i] = Qw[i].w;
    gaussian_elimination(A, x, b);
    for ( size_t i = 0; i < x.size(); ++i ) {
	Pw[i].w = x[i];
    }

    return SUCCESS;
}
    

//****f*
NurbsSurface skinned_surface(const vector<Nurbs> &C, int q)
{
/** Given a sequence of section curves, form a skinned surface.
 *
 * This function implements the skinning procedure described
 * in Piegl and Tiller (1997) for rational section curves.
 * See Section 10.3 of their text.
 *
 * This function assumes that all the curves are of the same
 * degree, but may have different knot vectors.  The algorithm
 * itself is not limited to curves of the same degree, but I haven't
 * implemented that here; I haven't had the need.
 *
 * Reference:
 * Piegl and Tiller (1997)
 * The NURBS Book, 2nd edition
 * Springer-Verlag, Berlin
 *
 * Input:
 * C    -- vector of section curves (assumed to be of same degree)
 * q    -- chosen degree for the interpolating v curves
 *
 * Return:
 * N    -- a NurbsSurface whicn "skins" the section curves
 *
 **/

    if ( q >= (int) C.size() ) {
	cout << "skinned_surface():\n";
	cout << "The chosen degree q= " << q << " is too large.\n";
	cout << "It should be less than the number of supplied curves: " << C.size() << endl;
	cout << "Bailing out!\n";
	exit(1);
    }

    // 1. Check that all curves have the same degree.
    //    cout << "1. Check that all curves have the same degree.\n";
    int p_ref = C[0].p;
    for ( size_t k = 1; k < C.size(); ++k ) {
	if ( C[k].p != p_ref ) {
	    cout << "skinned_surface():\n";
	    cout << "Error, not all supplied section curves share the same degree.\n";
	    cout << "deg(C[0])= " << p_ref << " deg(C[" << k << "])= " << C[k].p << endl;
	    cout << "Bailing out!\n";
	    exit(1);
	}
    }

    // 2. Make all curves with a common knot vector.
    //    cout << "2. Make all curves with a common knot vector.\n";
    vector<double> U(C[0].U);
    for ( size_t k = 1; k < C.size(); ++k ) {
	size_t i = 0;
	while ( i < C[k].U.size() ) {
	    double ui = C[k].U[i];
	    size_t count_ui = count(C[k].U.begin(), C[k].U.end(), ui);

	    vector<double>::iterator found;
	    found = find(U.begin(), U.end(), ui);

	    if ( found == U.end() ) {
		// Easy case: ui is not in U. So add count_ui lots of
		// ui to U.
		for ( size_t j = 0; j < count_ui; ++j ) {
		    U.push_back(ui);
		}
	    }
	    else {
		// Only add if multiplicity is greater.
		size_t count_ui2 = count(U.begin(), U.end(), ui);
		int diff = count_ui - count_ui2;
		if ( diff > 0 ) {
		    for ( int j = 0; j < diff; ++j ) {
			U.push_back(ui);
		    }
		}
	    }
	    i = i + count_ui;
	}
    }

    sort(U.begin(), U.end());

    // Knot refine based on difference in sets of knots
    vector<double> X;
    vector<Nurbs> C1;
    
    for ( size_t k = 0; k < C.size(); ++k ) {
	X.clear();
	for ( size_t i = 0; i < U.size(); ++i ) {
	    vector<double>::const_iterator found;
	    found = find(C[k].U.begin(), C[k].U.end(), U[i]);
	    if ( found == C[k].U.end() ) {
		// U[j] not found, so insert.
		X.push_back(U[i]);
	    }
	}
	C1.push_back(C[k].knot_refinement(X));
    }


    // 3. Setup {v_k}
    //    cout << "3. Setup {v_k}\n";

    int n = C1[0].P.size()-1;
    vector<double> vk(C1.size(), 0.0);
    vector<double> d(n, 0.0);

    for ( int i = 0; i <= n; ++i ) {
	for ( size_t k = 1; k < C1.size(); ++k ) {
	    d[i] += abs(C1[k].Pw[i] - C1[k-1].Pw[i]);
	    //d[i] += vabs(C1[k].P[i] - C1[k-1].P[i]);
	}
    }

    vk[0] = 0.0;
    vk[C1.size()-1] = 1.0;

    for ( size_t k = 1; k < C1.size()-1; ++k ) {
	double sum = 0.0;
	
	for ( int i = 0; i <= n; ++i ) {
	    sum += abs(C1[k].Pw[i] - C1[k-1].Pw[i])/d[i];
	    //sum += vabs(C1[k].P[i] - C1[k-1].P[i])/d[i];
	}

	vk[k] = vk[k-1] + (1.0/(n + 1.0))*sum;
    }

    // 4. Knot vector
    //    Using Eq. 9.8
    //    cout << "4. Knot vector\n";
    n = C1.size() - 1;
    int m = q + n + 1;

    vector<double> V(m+1, 0.0);

    for ( int i = 0; i <= q; ++i ) {
	V[i] = 0.0;
    }

    for ( int i = m-q; i <= m; ++i ) {
	V[i] = 1.0;
    }

    for ( int j = 1; j <= n - q; ++j ) {
	V[j+q] = 0.0;
	for ( int i = j; i <= j+q-1; ++i ) {
	    V[j+q] += vk[i];
	}
	V[j+q] *= 1.0/q;
    }

    // 5. Now ready to perform n+1 interpolations
    //    cout << "5. Now ready to perform n+1 interpolations\n";
    n = C1[0].P.size()-1;
    vector<vector<Mapped_point> > Pw(n+1);
    vector<Mapped_point> Pwi;
    vector<Mapped_point> Qw;

    for ( int i = 0; i <= n; ++i ) {

	Pw[i].resize(C1.size());
	Qw.clear();
	for ( size_t k = 0; k < C1.size(); ++k ) {
	    Qw.push_back(C1[k].Pw[i]);
	}

	interp_homogeneous_points(Qw, q, vk, V, Pwi);

	for ( size_t k = 0; k < C1.size(); ++k ) {
	    Pw[i][k] = Pwi[k];
	}
    }

    return NurbsSurface(Pw, p_ref, U, q, V);

}
