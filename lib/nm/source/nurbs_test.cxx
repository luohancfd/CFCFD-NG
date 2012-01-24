// Author: Rowan J. Gollan
// Date: 10-Sep-2008

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include "nurbs.hh"

using namespace std;

int test_find_span(ofstream &out)
/** Tests the implementation of find_span() in nurbs.cxx.
 *
 *  This is part of Ex2.3 from Piegl and Tiller (1997).
 **/
{
    int p = 2;
    vector<double> U(11, 0);
    U[0] = 0; U[1] = 0; U[2] = 0; U[3] = 1;
    U[4] = 2; U[5] = 3; U[6] = 4; U[7] = 4;
    U[8] = 5; U[9] = 5; U[10] = 5;
    double u = 5.0/2.0;
    int m = U.size() - 1; // m+1 == no. of knots.
    int n = m - p - 1;

    int i = find_span(u, n, p, U);

    bool result = ( i == 4 );

    out << "Test    'nurbs.cxx: find_span'" << endl;
    out << "Type    'function'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;

    return 0;
}

int test_basis_funs(ofstream &out)
/** Tests the implementation of basis_funs in nurbs.cxx.
 *
 *  This is the rest of Ex2.3 from Piegl and Tiller (1997).
 **/
{
    const double tol = 1.0e-9;
    
    int p = 2;
    vector<double> U(11, 0);
    U[0] = 0; U[1] = 0; U[2] = 0; U[3] = 1;
    U[4] = 2; U[5] = 3; U[6] = 4; U[7] = 4;
    U[8] = 5; U[9] = 5; U[10] = 5;
    double u = 5.0/2.0;
    int m = U.size() - 1; // m+1 == no. of knots.
    int n = m - p - 1;

    vector<double> N(p+1, 0.0);

    int i = 4;
    basis_funs(u, i, p, U, N);
    bool result = ( fabs(N[0] - 1./8.) <= tol &&
		    fabs(N[1] - 6./8.) <= tol &&
		    fabs(N[2] - 1./8.) <= tol );
    
    p = 1;
    N.resize(p+1);
    basis_funs(u, i, p, U, N);
    result = ( result &&
	       fabs(N[0] - 1./2.) <= tol &&
	       fabs(N[1] - 1./2.) <= tol );

    p = 0;
    N.resize(p+1);
    basis_funs(u, i, p, U, N);
    result = ( result &&
	       fabs(N[0] - 1.0) <= tol );

    
    out << "Test    'nurbs.cxx: basis_funs'" << endl;
    out << "Type    'function'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;

    return 0;
}

int test_nurbs_curve_point(ofstream &out)
/** Tests the implementation of nurbs_curve_point in nurbs.cxx.
 *
 *  This is Ex4.1 in Piegl and Tiller (1997).
 **/
{
    const double tol = 1.0e-9;
    int p = 2;

    vector<double> U(8, 0);
    U[0] = 0; U[1] = 0; U[2] = 0; U[3] = 1;
    U[4] = 2; U[5] = 3; U[6] = 3; U[7] = 3;
    
    vector<mapped_point> Pw;
    Pw.resize(5);
    Pw[0].wx = 0.0; Pw[0].wy = 0.0; Pw[0].wz = 0.0; Pw[0].w = 1.0;
    Pw[1].wx = 4.0; Pw[1].wy = 4.0; Pw[1].wz = 0.0; Pw[1].w = 4.0;
    Pw[2].wx = 3.0; Pw[2].wy = 2.0; Pw[2].wz = 0.0; Pw[2].w = 1.0;
    Pw[3].wx = 4.0; Pw[3].wy = 1.0; Pw[3].wz = 0.0; Pw[3].w = 1.0;
    Pw[4].wx = 5.0; Pw[4].wy = -1.0; Pw[4].wz = 0.0; Pw[4].w = 1.0;

    vector<double> C(3, 0.0);
    double u = 1.0;
    
    nurbs_curve_point(u, p, U, Pw, C);

    bool result = ( fabs(C[0] - 7./5.) <= tol &&
		    fabs(C[1] - 6./5.) <= tol );

    out << "Test    'nurbs.cxx: nurbs_curve_point'" << endl;
    out << "Type    'function'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;

    return 0;
}



int main()
{
    ofstream out;
    out.open("nurbs_test.result");
    if ( out.fail() ) {
	cout << "Error opening file nurbs_test.result\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    test_find_span(out);
    test_basis_funs(out);
    test_nurbs_curve_point(out);

    out.close();

    return 0;
}
