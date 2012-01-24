/* \file nelmin.hh
 * \ingroup nm
 * \brief Nelder-Mead simplex minimization of a nonlinear (multivariate) function.
 * \author PJ
 * \version 29-Jan-06 adapted from nelmin.py
 * 
 * The Python code had been adapted from the C-coded nelmin.c which was
 * adapted from the Fortran-coded nelmin.f which was, in turn, adapted
 * from the papers
 *
 *     J.A. Nelder and R. Mead (1965)
 *     A simplex method for function minimization.
 *     Computer Journal, Volume 7, pp 308-313.
 *
 *     R. O'Neill (1971)
 *     Algorithm AS47. Function minimization using a simplex algorithm.
 *     Applied Statistics, Volume 20, pp 338-345.
 *
 * and some examples are in
 *
 *    D.M. Olsson and L.S. Nelson (1975)
 *    The Nelder-Mead Simplex procedure for function minimization.
 *    Technometrics, Volume 17 No. 1, pp 45-51.
 *   
 * For a fairly recent and popular incarnation of this minimizer,
 * see the amoeba function in the famous "Numerical Recipes" text.
 * The programming interface is via the minimize() function; see below.
 */

#ifndef NELMIN_HH
#define NELMIN_HH

#include <vector> 
using namespace std;
#include "fobject.hh"

//-----------------------------------------------------------------------
// The public face of the minimizer...

// typedef double (*objectiveFn)(vector<double> &);
// Have replaced the pointer-to function with a function object.


/// \brief Locate a minimum of the objective function, f.
///
/// Input:
/// f       : user-supplied MultivariateFunction object f(x)
/// x       : vector of N coordinates
/// dx      : vector of N increments to apply to x when forming
///           the initial simplex.  Their magnitudes determine the size
///           and shape of the initial simplex.
/// tol     : the terminating limit for the standard-deviation
///           of the simplex function values.
/// maxfe   : maximum number of function evaluations that we will allow
/// n_check : number of steps between convergence checks
/// delta   : magnitude of the perturbations for checking a local minimum
///           and for the scale reduction when restarting
/// Kreflect, Kextend, Kcontract: coefficients for locating the new vertex
///
/// Output:
/// Returns a flag to indicate if convergence was achieved.
/// On return x contains the coordinates for the best x location,
/// corresponding to min(f(x)),
/// f_min : the function value at that point,
/// n_fe : the number of function evaluations and
/// n_restart : the number of restarts (with scale reduction).
int minimize( MultivariateFunction *f, 
	      vector<double> &x, 
	      double *f_min,
	      int *n_fe,
	      int *n_restart,
	      vector<double> *dx=0,
	      double tol=1.0e-6,
	      int max_fe=300, 
	      int n_check=20, 
	      double delta=0.001,
	      double Kreflect=1.0, 
	      double Kextend=2.0, 
	      double Kcontract=0.5);

//-----------------------------------------------------------------------
// Use a class to keep the data tidy and conveniently accessible...

/// \brief Stores the (nonlinear) simplex as a list of lists.
///
/// In an N-dimensional problem, each vertex is a list of N coordinates
/// and the simplex consists of N+1 vertices.
class NMSimplex {
    public:
    int N;
    vector<double> dx;
    vector<vector<double>*> vertex_list;
    vector<double> f_list;
    MultivariateFunction *f;
    int n_fe;
    int n_restarts;
    double Kreflect;
    double Kextend;
    double Kcontract;
    /// Set up the vertices about the user-specified vertex, x, and the step-sizes dx.
    NMSimplex( vector<double> &_x, vector<double> &_dx, MultivariateFunction *_f, 
	       double _Kreflect, double _Kextend, double _Kcontract);
    ~NMSimplex();
    void take_a_step();
    void rescale( double ratio );
    int lowest( int exclude=-1 );
    int highest( int exclude=-1 );
    double std_dev();
    vector<double> centroid( int exclude=-1 );
    void contract_about_one_point( int i_con );
    bool test_for_minimum( int i_min, double delta );
    void replace_vertex( int i, vector<double>* p, double fp );
};

// Create a new N-dimensional point as a weighting of points p1 and p2.
vector<double>* create_new_point( double c1, vector<double> &p1, 
				  double c2, vector<double> &p2);

#endif
