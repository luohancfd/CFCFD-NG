/* \file nelmin_test.cxx
 * \ingroup nm
 * \brief Try out the Nelder-Mead Simplex optimizer.
 * \author PJ
 * \version 29-Jan-06 adapted from nelmin.py
 */

#include <math.h>
#include <iostream>
#include <vector>
using namespace std;
#include "fobject.hh"
#include "nelmin.hh"

// Test objective function 1.
// x is expected to be a list of ccordinates.
// Returns a single float value.
class TestFunction1: public MultivariateFunction { 
public:
    double eval( vector<double> &x )
    {
	size_t n = x.size();
	double sum = 0.0;
	for ( size_t i = 0; i < n; ++i ) {
	    sum += (x[i] - 1.0) * (x[i] - 1.0);
	}
	return sum;
    }
};
 
// Test objective function 2.
// Example 3.3 from Olsson and Nelson.
class TestFunction2: public MultivariateFunction { 
public:
    double eval( vector<double> &x )
    {
	double x1 = x[0]; double x2 = x[1];   // rename to match the paper
	if ( (x1 * x1 + x2 * x2) > 1.0 ) {
	    return 1.0e38;
	} else {
	    double yp = 53.69 + 7.26 * x1 - 10.33 * x2 + 7.22 * x1 * x1 
		+ 6.43 * x2 * x2 + 11.36 * x1 * x2;
	    double ys = 82.17 - 1.01 * x1 - 8.61 * x2 + 1.40 * x1 * x1 
		- 8.76 * x2 * x2 - 7.20 * x1 * x2;
	    return -yp + fabs(ys - 87.8);
	}
    }
};
 
// Test objective function 3.
// Example 3.5 from Olsson and Nelson; least-squares.
class TestFunction3: public MultivariateFunction { 
public:
    double eval( vector<double> &z )
    {
	double x[] = {0.25, 0.50, 1.00, 1.70, 2.00, 4.00};
	double y[] = {0.25, 0.40, 0.60, 0.58, 0.54, 0.27};
	double a1 = z[0]; double a2 = z[1]; 
	double alpha1 = z[2]; double alpha2 = z[3];
	double sum_residuals = 0.0;
	for ( size_t i = 0; i < 6; ++i ) {
	    double t = x[i];
	    double eta = a1 * exp(alpha1 * t) + a2 * exp(alpha2 * t);
	    double r = y[i] - eta;
	    sum_residuals += r * r;
	}
	return sum_residuals;
    }
};
   
// -------------------------------------------------------------------

int main() {
    cout << "Begin nelmin self-test..." << endl;

    cout << "---------------------------------------------------" << endl;
    cout << "test 1: simple quadratic with zero at (1,1,...)" << endl;
    vector<double> x = vector<double>(3);
    double fx;
    int nfe, nres;
    TestFunction1 test_fun_1 = TestFunction1();
    bool conv_flag = minimize( &test_fun_1, x, &fx, &nfe, &nres);
    cout << "x= (" << x[0] << " " << x[1] << " " << x[2] << ")" << endl;
    cout << "fx=" << fx << endl;
    cout << "convergence-flag=" << conv_flag << endl;
    cout << "number-of-fn-evaluations=" << nfe << endl;
    cout << "number-of-restarts=" << nres << endl;

    cout << "---------------------------------------------------" << endl;
    cout << "test 2: Example 3.3 in Olsson and Nelson f(0.811,-0.585)=-67.1" << endl;
    TestFunction2 test_fun_2 = TestFunction2();
    vector<double> x2 = vector<double>(2);
    vector<double> dx2 = vector<double>(2, 0.5);
    conv_flag = minimize( &test_fun_2, x2, &fx, &nfe, &nres, &dx2, 1.0e-4);
    cout << "x= (" << x2[0] << " " << x2[1] << ")" << endl;
    cout << "fx=" << fx << endl;
    cout << "convergence-flag=" << conv_flag << endl;
    cout << "number-of-fn-evaluations=" << nfe << endl;
    cout << "number-of-restarts=" << nres << endl;

    cout << "---------------------------------------------------" << endl;
    cout << "test 3: Example 3.5 in Olsson and Nelson, nonlinear least-squares" << endl;
    cout << "f(1.801, -1.842, -0.463, -1.205)=0.0009" << endl;
    TestFunction3 test_fun_3 = TestFunction3();
    vector<double> x3 = vector<double>(4, 1.0);
    x3[2] = -0.5; x3[3] = -2.5;
    vector<double> dx3 = vector<double>(4, 0.1);
    conv_flag = minimize( &test_fun_3, x3, &fx, &nfe, &nres, &dx3, 1.0e-9, 800);
    cout << "x= (" << x3[0] << " " << x3[1] << " " << x3[2] << " " << x3[3] << ")" << endl;
    cout << "fx=" << fx << endl;
    cout << "convergence-flag=" << conv_flag << endl;
    cout << "number-of-fn-evaluations=" << nfe << endl;
    cout << "number-of-restarts=" << nres << endl;

    cout << "---------------------------------------------------" << endl;
    cout << "Done." << endl;
    return 0;
}
