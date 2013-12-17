/*  \file ode_test.cxx
 *  \brief A testing program for the ODE suite
 *
 *  \author Rowan J Gollan
 *  \version 21-Feb-2006
 **/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <vector>
#include "no_fuss_linear_algebra.hh"
#include "ode_system.hh"
#include "ode_solver.hh"
#include "ode_step.hh"
using namespace std;

class Sys1 : public OdeSystem {
public:
    Sys1( int ndim, bool system_test );
    Sys1( const Sys1 &s );
    virtual ~Sys1();

    int eval( const vector<double> &y, vector<double> &ydot );
};

Sys1::Sys1( int ndim, bool system_test )
    : OdeSystem( ndim, system_test ) {}

Sys1::Sys1( const Sys1 &s )
    : OdeSystem( s.ndim_, s.apply_system_test_ ) {}

Sys1::~Sys1() {}

int Sys1::eval( const vector<double> &y, vector<double> &ydot )
{
    ydot[0] = -0.5 * y[0];
    ydot[1] = y[0] + -5.0 * y[1];

    return 0;
}

class Sys2 : public OdeSystem {
public:
    Sys2( int ndim, bool system_test );
    Sys2( const Sys2 &s );
    virtual ~Sys2();

    int eval( const vector<double> &y, vector<double> &ydot );
    int eval_split( const vector<double> &y, vector<double> &p,
		    vector<double> &q );
};

Sys2::Sys2( int ndim, bool system_test )
    : OdeSystem( ndim, system_test ) {}

Sys2::Sys2( const Sys2 &s )
    : OdeSystem( s.ndim_, s.apply_system_test_ ) {}

Sys2::~Sys2() {}

int Sys2::eval( const vector<double> &y, vector<double> &ydot )
{
    ydot[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    ydot[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * pow(y[1], 2);
    ydot[2] = 3.0e7 * pow(y[1], 2);

    return 0;
}

int Sys2::eval_split( const vector<double> &y, vector<double> &q,
		      vector<double> &L )

{
    q[0] = 1.0e4 * y[1] * y[2];
    L[0] = 0.04 * y[0];

    q[1] = 0.04 * y[0];
    L[1] = 1.0e4 * y[1] * y[2] + 3.0e7 * pow(y[1], 2);

    q[2] = 3.0e7 * pow(y[1], 2);
    L[2] = 0.0;

    return 0;
}


void printUsage()
{
    cout << "Usage: ode_test.x [--verbose|--fp-only|--data]" << endl;
    cout << "Options: \n";
    cout << "    --help      : Print this help.\n";
    cout << "    --verbose   : Print tests in human-readable form (default action)\n";
    cout << "    --fp-only   : Same as verbose but only floating-point values\n";
    cout << "                  are printed.  This is useful for regression comparisons.\n";
    cout << "    --data      : Generate data over large T range for plotting in an output file\n";
    return;
}

void run_standard_test( bool verbose );
void generate_data();

int main( int argc, char *argv[] )
{
    if( argc > 2 ) {
	printUsage();
	return 0;
    }

    if( argc == 1 ) {
	run_standard_test( true );
	return 0;
    }

    string action( argv[1] );

    if( action == "--help" ) {
	printUsage();
	return 0;
    }
    else if( action == "--verbose" ) {
	run_standard_test( true );
	return 0;
    }
    else if( action == "--fp-only" ) {
	run_standard_test( false );
	return 0;
    }
    else if( action == "--data" ) {
	generate_data();
	return 0;
    }
    else {
	cout << "Uknown option: " << action << endl;
	printUsage();
	return 0;
    }

    return 0;

}

void run_standard_test( bool verbose )
{

    cout << setprecision(10) << showpoint;

    if( verbose ) {
	cout << "Begin test of OdeSolver, OdeSystem, OdeStep classes...\n";

	cout << "------------------------------------------------------\n";
	cout << " Test 1: A 2x2 ODE system y' = Ay   \n";
	cout << "------------------------------------------------------\n";
	cout << "A = | -0.5   0.0 |  y(0) = | 1.0 |                    \n";
	cout << "    |  1.0  -5.0 |         | 0.0 |                  \n\n";
	cout << "Using a stepsize= 0.01, we solve for y(1.0).          \n";
	cout << endl;
	cout << "Analytical result:                                    \n";
	cout << "y(1.0) = | 0.6065306 |                                 \n";
	cout << "         | 0.1332873 |                               \n\n";
    }
  
    Sys1 system1( 2, false );
    OdeSolver test1( string("test1"), 2, string("euler"), 4, 5.0 );
    vector<double> yin(2);
    vector<double> yout(2);

    yin[0] = 1.0; yin[1] = 0.0;
    copy_vector(yin, yout);
    double h = 0.01;

    test1.solve_over_interval( system1, 0.0, 1.0, &h, yin, yout );

    if( verbose ) {
	cout << "Euler method:\n"
	     << "y(1.0) = | " << yout[0] << " | " << endl
	     << "         | " << yout[1] << " |\n\n";
    }
    else {
	cout << yout[0] << "  " << yout[1] << endl;
    }
    

    ModEulerStep mod_step("mod_step", 2);
    test1.set_step( &mod_step );
    test1.solve_over_interval( system1, 0.0, 1.0, &h, yin, yout );

    if( verbose ) {
	cout << "Modified Euler method:\n"
	     << "y(1.0) = | " << yout[0] << " | " << endl
	     << "         | " << yout[1] << " |\n\n";
    }
    else {
	cout << yout[0] << "  " << yout[1] << endl;
    }

    RKFStep rkf_step("rkf_step", 2, 1.0e-8);
    test1.set_step( &rkf_step );
    
    test1.solve_over_interval( system1, 0.0, 1.0, &h, yin, yout );

    if( verbose ) {
	cout << "Runge-Kutta-Fehlberg method:\n"
	     << "y(1.0) = | " << yout[0] << " | " << endl
	     << "         | " << yout[1] << " |\n\n";
	cout << "Timestep at end: " << h << endl;
    }
    else {
	cout << yout[0] << "  " << yout[1] << endl;
    }
    cout << endl;

    h = 0.01;

    DP853Step dp853_step("dp853_step", 2, 1.0e-8);
    test1.set_step( &dp853_step );
    
    test1.solve_over_interval( system1, 0.0, 1.0, &h, yin, yout );

    if( verbose ) {
	cout << "Dormand-Prince 8(53) RK method:\n"
	     << "y(1.0) = | " << yout[0] << " | " << endl
	     << "         | " << yout[1] << " |\n\n";
	cout << "Timestep at end: " << h << endl;
    }
    else {
	cout << yout[0] << "  " << yout[1] << endl;
    }
    cout << endl;

    h = 0.01;

    QssStep qss_step("qss_step", 2, 4, 1.0e-5, 2.0);
    test1.set_step( &qss_step );
    
    test1.solve_over_interval( system1, 0.0, 1.0, &h, yin, yout );

    if( verbose ) {
	cout << "alpha-QSS method:\n"
	     << "y(1.0) = | " << yout[0] << " | " << endl
	     << "         | " << yout[1] << " |\n\n";
	cout << "Timestep at end: " << h << endl;
    }
    else {
	cout << yout[0] << "  " << yout[1] << endl;
    }

    
    if( verbose ) {
	cout << "------------------------------------------------------\n";
	cout << " Test 2: CVODE example (stiff problem)  \n";
	cout << "------------------------------------------------------\n";
	cout << "Solve the following chemical kinetics problem..       \n";
	cout << "dy1/dt = -0.04 * y1 + 1.0e4*y2*y3\n";
	cout << "dy2/dt = 0.04 * y1 - 1.e4*y2*y3 - 3.0e7*(y2)^2\n";
	cout << "dy3/dt = 3.0e7*(y2)^2\n\n";
	cout << "on the interval t = 0.0 to t = 40.0 with\n\n";
	cout << "y0 = | 1.0 |\n";
	cout << "     | 0.0 |\n";
	cout << "     | 0.0 |\n\n";
	cout << "Table of results:\n";
	cout << "========================================================\n";
	cout << "     t     |     y[0]     |     y[1]     |      y[2]    \n";
	cout << "========================================================\n";
	cout << " CVODE reference answer                                  \n";
	cout << "    0.4    | 9.851641e-01 | 3.386242e-05 | 1.480205e-02 \n";
	cout << "    4.0    | 9.055097e-01 | 2.240338e-05 | 9.446793e-02 \n";
	cout << "   40.0    | 7.158016e-01 | 9.185045e-06 | 2.841893e-01 \n";
    }
  
    Sys2 system2( 3, false );
    OdeSolver test2( string("test2"), 3, string("euler"), 4,
		     4.0, 0.001, 0.33);
    yin.resize(3);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    yout.resize(3);
    h = 0.0002;

    test2.solve_over_interval( system2, 0.0, 4.0e-1, &h, yin, yout );

    if( verbose ) {

	cout << "\n Euler method (probably not great for a stiff problem)  \n";
	printf( "    0.4    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 0.4, 4.0, &h, yin, yout );

    if( verbose ) {

	printf( "    4.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }
    
    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 4.0, 40.0, &h, yin, yout );

    if( verbose ) {

	printf( "   40.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    ModEulerStep mod_step2("mod_step", 3);
    test2.set_step( &mod_step2 );

    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    
    test2.solve_over_interval( system2, 0.0, 4.0e-1, &h, yin, yout );

    if( verbose ) {

	cout << "\n Modified Euler method \n";
	printf( "    0.4    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 0.4, 4.0, &h, yin, yout );

    if( verbose ) {

	printf( "    4.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 4.0, 40.0, &h, yin, yout );

    if( verbose ) {

	printf( "   40.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    RKFStep rkf_step2("rkf_step", 3, 1.0e-12);
    test2.set_step( &rkf_step2 );

    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    
    test2.solve_over_interval( system2, 0.0, 4.0e-1, &h, yin, yout );

    if( verbose ) {

	cout << "\n Runge-Kutta-Fehlberg method \n";
	printf( "    0.4    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 0.4, 4.0, &h, yin, yout );

    if( verbose ) {

	printf( "    4.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 4.0, 40.0, &h, yin, yout );

    if( verbose ) {

	printf( "   40.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);
	cout << "     --- final  timestep= " << h << " ---" << endl;

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }


    h = 0.0002;

    DP853Step dp853_step2("dp853_step", 3, 1.0e-12);
    test2.set_step( &dp853_step2 );

    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    
    test2.solve_over_interval( system2, 0.0, 4.0e-1, &h, yin, yout );

    if( verbose ) {

	cout << "\n Dormand-Prince 8(53) RK method \n";
	printf( "    0.4    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 0.4, 4.0, &h, yin, yout );

    if( verbose ) {

	printf( "    4.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 4.0, 40.0, &h, yin, yout );

    if( verbose ) {

	printf( "   40.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);
	cout << "     --- final  timestep= " << h << " ---" << endl;

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    h = 0.0002;

    QssStep qss_step2("qss_step", 3, 10, 1.0e-6, 5.0);
    test2.set_step( &qss_step2 );

    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    
    test2.solve_over_interval( system2, 0.0, 4.0e-1, &h, yin, yout );

    if( verbose ) {

	cout << "\n alpha-QSS method \n";
	printf( "    0.4    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 0.4, 4.0, &h, yin, yout );

    if( verbose ) {

	printf( "    4.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }

    copy_vector(yout, yin);
    test2.solve_over_interval( system2, 4.0, 40.0, &h, yin, yout );

    if( verbose ) {

	printf( "   40.0    | %12.6e | %12.6e | %12.6e \n", yout[0], yout[1], yout[2]);
	cout << "     --- final  timestep= " << h << " ---" << endl;

    }
    else {
	cout << yout[0] << "  " << yout[1] << "  " << yout[2] << endl;
    }


    return;
}

    
void generate_data()
{
    return;
}

    
