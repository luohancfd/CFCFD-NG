// ridder_test.cxx
// Solve a nonlinear equation f(x)=0 using the method of Ridder.
// Peter J. 
// mech3750 demo code 12-Mar-2014

#include <iostream>
#include <functional>
#include <math.h>
#include "ridder.hh"

using namespace std;

double test_fun_1(double x) {
    return pow(x,3) + pow(x,2) - 3*x - 3;
}

double test_fun_2(double x, double a) {
    return a*x + sin(x) - exp(x);
}

int main() {
    cout << "Begin self-test of Ridder's method..." << endl;
    cout << endl;
    cout << "Example from Gerald and Wheatley, p. 45" << endl;
    cout << "Solve f(x) = x^3 + x^2 - 3x -3 = 0 with initial" << endl;
    cout << "guesses of x0 = 1 and x1 = 2." << endl;
    cout << "Final result x = " << solve(test_fun_1, 1, 2) << endl;
    cout << "Gerald and Wheatley report x = 1.732051" << endl;
    cout << endl;
    //
    cout << "Example from Gerald and Wheatley, p.45 also, " 
	 << "but using lambda expression." << endl;
    cout << "Solve f(x) = 3*x + sin(x) - e^x = 0 with initial" << endl;
    cout << "guesses of x0 = 0 and x1 = 1." << endl;
    double my_a = 3.0;
    auto test_fun_3 = [my_a] (double x) -> double
	{ return test_fun_2(x,my_a); }; 
    cout << "Final result x = " 
	 <<  solve(test_fun_3, 0, 1) << endl;
    cout << "Gerald and Wheatley report x = 0.3604217" << endl;
    cout << endl;
    //
    cout << "Bracket a root of the second example." << endl;
    double x1 = 0.4;
    double x2 = 0.5;
    int result_flag = bracket(test_fun_3, x1, x2);
    cout << "result_flag=" << result_flag 
	 << " x1=" << x1 
	 << " x2=" << x2 << endl;
    cout << "Done." << endl;
    //
    return 0;
}
