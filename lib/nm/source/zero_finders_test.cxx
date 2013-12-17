/** \file zero_systems_test.cxx
 * \brief Testing of the zero_finders module
 *
 * \author Rowan J Gollan
 * \version 24-Apr-2006
 *
 **/

#include <cmath>
#include <iostream>
#include <vector>

#include "no_fuss_linear_algebra.hh"
#include "zero_system.hh"
#include "zero_finders.hh"

using namespace std;

class SystemOne : public ZeroSystem {
public:
    SystemOne();
    ~SystemOne();

    int f( const vector<double> &y, vector<double> &G );
    int Jac( const vector<double> &y, Valmatrix &dGdy );

};

SystemOne::SystemOne()
    : ZeroSystem() {}

SystemOne::~SystemOne() {}

int SystemOne::f( const vector<double> &y, vector<double> &G )
{
    G[0] = 4.0 - y[0]*y[0] - y[1]*y[1];
    G[1] = 1.0 - exp(y[0]) - y[1];

    return 0;
}

int SystemOne::Jac( const vector<double> &y, Valmatrix &dGdy )
{
    dGdy.set(0, 0, -2.0 * y[0]);  dGdy.set(0, 1, -2.0 * y[1]);
    dGdy.set(1, 0, -1.0 * exp(y[0])); dGdy.set(1, 1, -1.0);

    return 0;

}


class SystemTwo : public ZeroSystem {
public:
    SystemTwo();
    ~SystemTwo();

    int f( const vector<double> &y, vector<double> &G );
    int Jac( const vector<double> &y, Valmatrix &dGdy );

};

SystemTwo::SystemTwo()
    : ZeroSystem() {}

SystemTwo::~SystemTwo() {}

int SystemTwo::f( const vector<double> &y, vector<double> &G )
{
    G[0] = exp(y[0]) - y[1];
    G[1] = y[0] * y[1] - exp(y[0]);

    return 0;
}

int SystemTwo::Jac( const vector<double> &y, Valmatrix &dGdy )
{
    dGdy.set(0,0, exp(y[0]));   dGdy.set(0,1, -1.0);
    dGdy.set(1,0, y[1] - exp(y[0])); dGdy.set(1,1, y[0]);
    
    return 0;
}


class SystemThree : public ZeroSystem {
public:
    SystemThree();
    ~SystemThree();

    int f( const vector<double> &y, vector<double> &G );
    int Jac( const vector<double> &y, Valmatrix &dGdy );

};

SystemThree::SystemThree()
    : ZeroSystem() {}

SystemThree::~SystemThree() {}

int SystemThree::f( const vector<double> &y, vector<double> &G )
{
    G[0] = y[0] * y[0] - 2.0 * y[0] - y[1] + 0.5;
    G[1] = y[0] * y[0] + 4.0 * y[1] * y[1] - 4.0;

    return 0;
}

int SystemThree::Jac( const vector<double> &y, Valmatrix &dGdy )
{

    dGdy.set(0,0, 2.0 * y[0] - 2.0);   dGdy.set(0,1, -1.0);
    dGdy.set(1,0, 2.0 * y[0]);         dGdy.set(1,1, 8.0 * y[1]);
    
    return 0;
}




int main()
{
    cout << "Begin zero_finders_test.x...\n\n";
    
    cout << "-----------------------------------------\n";
    cout << "--- Test case 1: 2x2 nonlinear system ---\n";
    cout << "-----------------------------------------\n";

    cout << "This is from page 175 of Gerald and Wheatley.\n";
    cout << "Solve for x and y the set of equations:\n\n";
    cout << "x^2 + y^2 = 4 \n";
    cout << "e^x + y = 1 \n";
    cout << "Take as an initial guess: x = 1.0, y = -1.7\n";
    NewtonRaphsonZF NRSolver( 2, 1.0e-6, 10 );
    SystemOne Sys1;
    vector<double> y_guess( 2 );
    y_guess[0] = 1.0; y_guess[1] = -1.7;

    vector<double> y_out( 2 );
    NRSolver.solve( Sys1, y_guess, y_out );

    cout << "After the call to the Newton-Raphson solver...\n";
    cout << "NRSolver.solve()\n";
    cout << "y_out= \n";
    print_vector(y_out);
   
    cout << "The answer given in Gerald and Wheatley is :\n";
    cout << "x = 1.004167, y = -1.729635\n";

    cout << "-----------------------------------------\n";
    cout << "--- Test case 2: 2x2 nonlinear system ---\n";
    cout << "-----------------------------------------\n";

    cout << "This is from page 177 of Gerald and Wheatley.\n";
    cout << "Solve for x and y the set of equations:\n\n";
    cout << "e^x - y   = 0 \n";
    cout << "x y - e^x = 0 \n";
    cout << "Take as an initial guess: x = 0.95, y = 2.7\n";
    SystemTwo Sys2;
    y_guess[0] = 0.95; y_guess[1] = 2.7;
    NRSolver.solve( Sys2, y_guess, y_out );

    cout << "After the call to the Newton-Raphson solver...\n";
    cout << "NRSolver.solve()\n";
    cout << "y_out= \n";
    print_vector(y_out);
   
    cout << "The answer given in Gerald and Wheatley is :\n";
    cout << "x = 1.000000, y = 2.718282\n";

    cout << "-----------------------------------------\n";
    cout << "--- Test case 3: 2x2 nonlinear system ---\n";
    cout << "-----------------------------------------\n";

    cout << "This is Example 2.20 from Mathews.\n";
    cout << "Solve for x and y the set of equations:\n\n";
    cout << "x^2 - 2x - y + 0.5 \n";
    cout << "x^2 + 4y^2 - 4 = 0 \n";
    cout << "Take as an initial guess: x = 2.00, y = 0.25\n";
    y_guess[0] = 2.00; y_guess[1] = 0.25;
    SystemThree Sys3;
    NRSolver.solve( Sys3, y_guess, y_out );

    cout << "After the call to the Newton-Raphson solver...\n";
    cout << "NRSolver.solve()\n";
    cout << "y_out= \n";
    print_vector(y_out);

    cout << "The answer given in Mathews is:\n";
    cout << "x = 1.900677, y = 0.311219\n";

    cout << "Done.\n";
    

    return 0;
}
    

