/* \file nelmin_test.cxx
 * \ingroup nm
 * \brief Implementation of Nelder-Mead simplex minimization.
 * \author PJ
 * \version 29-Jan-06 adapted from nelmin.py
 */

#include <math.h>
#include <iostream>
#include <vector>
using namespace std;
#include "nelmin.hh"

//-----------------------------------------------------------------------
// The public face of the minimizer...

int minimize( MultivariateFunction *f, vector<double> &x, 
	      double *f_min, int *n_fe, int *n_restart,
	      vector<double> *dx, double tol, int max_fe, 
	      int n_check, double delta,
	      double Kreflect, double Kextend, double Kcontract)
{
    bool converged = false;
    int N = x.size();
    bool did_allocate_dx = false;
    if ( dx == 0 ) {
	dx = new vector<double>(N, 0.1);
	did_allocate_dx = true;
    }
    NMSimplex smplx = NMSimplex(x, *dx, f, Kreflect, Kextend, Kcontract);
    *n_restart = 0;

    while ( !converged && smplx.n_fe < max_fe ) {
	// Take some steps and then check for convergence.
        for ( int i = 0; i < n_check; ++i ) {
            smplx.take_a_step();
	    // Pick out the current best vertex.
	    int i_best = smplx.lowest();
	    for ( int j = 0; j < N; ++j ) x[j] = (*smplx.vertex_list[i_best])[j];
	    *f_min = smplx.f_list[i_best];
	    // Check the scatter of vertex values to see if we are
	    // close enough to call it quits.
	    if ( smplx.std_dev() < tol ) {
		// All of the points are close together but we need to
		// test more carefully to see if we are at a true minimum.
		converged = smplx.test_for_minimum(i_best, delta);
		if ( !converged ) {
		    // The function evaluations are all very close together
		    // but we are not at a true minimum; rescale the simplex.
		    smplx.rescale(delta);
		    *n_restart += 1;
		}
	    }
	}
    }
    *n_fe = smplx.n_fe;
    if ( did_allocate_dx ) delete dx;
    return converged;
}


//-----------------------------------------------------------------------
// Use a class to keep the data tidy and conveniently accessible...

NMSimplex::NMSimplex( vector<double> &_x, vector<double> &_dx, MultivariateFunction *_f,
		      double _Kreflect, double _Kextend, double _Kcontract )
    : N(_x.size()), dx(_dx), 
      vertex_list(vector<vector<double>*>()), f_list(vector<double>()), 
      f(_f), n_fe(0),
      Kreflect(_Kreflect), Kextend(_Kextend), Kcontract(_Kcontract)
{
    // Set up the vertices about the user-specified vertex, x,
    // and the set of step-sizes dx.
    // f is a user-specified objective function f(x).
    for (int i = 0; i <= N; ++i ) {
	vector<double> *p = new vector<double>(_x);
	if ( i >= 1 ) (*p)[i-1] += dx[i-1];
       	vertex_list.push_back(p);
	f_list.push_back(f->eval(*p));
	n_fe += 1;
    }
}
NMSimplex::~NMSimplex()
{
    for (int i = 0; i <= N; ++i ) {
	delete vertex_list[i];
    }
}

//------------------------------------------------------------------------
// The core of the minimizer...
void NMSimplex::take_a_step()
{
    // Try to move away from the worst point in the simplex.
    // The new point will be inserted into the simplex (in place).
    int i_low = lowest();
    int i_high = highest();
    vector<double>* x_high = vertex_list[i_high];
    double f_high = f_list[i_high];
    // Centroid of simplex excluding worst point.
    vector<double> x_mid = centroid(i_high);
    double f_mid = f->eval(x_mid);
    n_fe += 1;

    // First, try moving away from worst point by
    // reflection through centroid
    vector<double>* x_refl = create_new_point(1.0+Kreflect, x_mid, -Kreflect, *x_high);
    double f_refl = f->eval(*x_refl);
    n_fe += 1;
    if ( f_refl < f_mid ) {
        // The reflection through the centroid is good,
        // try to extend in the same direction.
        vector<double>* x_ext = create_new_point(Kextend, *x_refl, 1.0-Kextend, x_mid);
        double f_ext = f->eval(*x_ext);
        n_fe += 1;
	if ( f_ext < f_refl ) {
            // Keep the extension because it's best.
            replace_vertex(i_high, x_ext, f_ext);
	    delete x_refl;  // clean up temporary data
	    return;
	} else {
            // Settle for the original reflection.
            replace_vertex(i_high, x_refl, f_refl);
	    delete x_ext;
	    return;
	}
    } else {
        // The reflection is not going in the right direction, it seems.
        // See how many vertices are better than the reflected point.
        int count = 0;
        for ( int i = 0; i <= N; ++i ) {
            if ( f_list[i] > f_refl )  count += 1;
	}
	if ( count <= 1 ) {
            // Not too many points are higher than the original reflection.
            // Try a contraction on the reflection-side of the centroid.
            vector<double>* x_con = create_new_point(1.0-Kcontract, x_mid, Kcontract, *x_high);
            double f_con = f->eval(*x_con);
            n_fe += 1;
	    if ( f_con < f_high ) {
                // At least we haven't gone uphill; accept.
                replace_vertex(i_high, x_con, f_con);
	    } else {
                // We have not been successful in taking a single step.
                // Contract the simplex about the current lowest point.
                contract_about_one_point(i_low);
		delete x_con;
	    }
	    delete x_refl; // We end up not using the original reflection.
	    return;
	} else {
            // Retain the original reflection because there are many
            // vertices with higher values of the objective function.
            replace_vertex(i_high, x_refl, f_refl);
	    return;
	}
    }
}

//-------------------------------------------------------------------------
// Utility functions.

// Pick out the current minimum and rebuild the simplex about that point.
void NMSimplex::rescale( double ratio )
{
    int i;
    int i_min = lowest();
    for ( i = 0; i < N; ++i ) dx[i] *= ratio;
    vector<double> x = vector<double>(*vertex_list[i_min]); // save to use below
    for ( i = 0; i <= N; ++i ) {
	for ( int j = 0; j < N; ++j ) (*vertex_list[i])[j] = x[j];
	if ( i >= 1 ) (*vertex_list[i])[i-1] += dx[i-1];
	f_list[i] = f->eval(*vertex_list[i]);
	n_fe += 1;
    }
    return;
}

// Returns the index of the lowest vertex, excluding the one specified.
int NMSimplex::lowest( int exclude ) 
{
    int indx;
    double lowest_f_value;
    if ( exclude == 0 ) indx = 1; else indx = 0;
    lowest_f_value = f_list[indx];
    for ( int i = 0; i <= N;  ++i ) {
	if ( i == exclude ) continue;
	if ( f_list[i] < lowest_f_value ) {
	    lowest_f_value = f_list[i];
	    indx = i;
	}
    }
    return indx;
}

// Returns the index of the highest vertex, excluding the one specified.
int NMSimplex::highest( int exclude ) 
{
    int indx;
    double highest_f_value;
    if ( exclude == 0 ) indx = 1; else indx = 0;
    highest_f_value = f_list[indx];
    for ( int i = 0; i <= N;  ++i ) {
	if ( i == exclude ) continue;
	if ( f_list[i] > highest_f_value ) {
	    highest_f_value = f_list[i];
	    indx = i;
	}
    }
    return indx;
}


// Returns the standard deviation of the vertex fn values.
double NMSimplex::std_dev() 
{
    int i;
    double sum = 0.0;
    for ( i = 0; i <= N; ++i ) {
	sum += f_list[i];
    }
    double mean = sum / (N + 1);
    sum = 0.0;
    for ( i = 0; i <= N; ++i ) {
	double diff = f_list[i] - mean;
	sum += diff * diff;
    }
    return sqrt(sum / N);
}

// Returns the centroid of all vertices excluding the one specified.
vector<double> NMSimplex::centroid( int exclude )
{
    int i, j;
    vector<double> xmid = vector<double>(N, 0.0);
    for ( i = 0; i <= N; ++i ) {
	if (i == exclude ) continue;
	for ( j = 0; j < N; ++j ) xmid[j] += (*vertex_list[i])[j];
    }
    for ( j = 0; j < N; ++j ) xmid[j] /= N;
    return xmid;
}

// Contract the simplex about the vertex i_con.    
void NMSimplex::contract_about_one_point( int i_con )
{
    vector<double> *p_con = vertex_list[i_con];
    for ( int i = 0; i <= N; ++i ) {
	if ( i == i_con ) continue;
	vector<double> *p = vertex_list[i];
	for ( int j = 0; j < N; ++j )
	    (*p)[j] = 0.5 * ((*p)[j] + (*p_con)[j]);
	f_list[i] = f->eval(*p);
	n_fe += 1;
    }
    return;
}

// Perturb the minimum vertex and check that it is a local minimum.
bool NMSimplex::test_for_minimum( int i_min, double delta )
{
    int is_minimum = 1;  // Assume it is true and test for failure.
    double f_min = f_list[i_min];
    double f_p;
    for ( int j = 0; j < N; ++j ) {
	// Check either side of the minimum, perturbing one coordinate at a time.
	vector<double> p = vector<double>( *vertex_list[i_min] );
	p[j] += dx[j] * delta;
	f_p = f->eval(p);
	n_fe += 1;
	if ( f_p < f_min ) { is_minimum = 0; break; }
	p[j] -= dx[j] * delta * 2;
	f_p = f->eval(p);
	n_fe += 1;
	if ( f_p < f_min ) { is_minimum = 0; break; }
    }
    return is_minimum;
}
   
// Create a new N-dimensional point as a weighting of points p1 and p2.
vector<double>* create_new_point( double c1, vector<double> &p1, 
				  double c2, vector<double> &p2)
{
    int n = p1.size();
    vector<double> *p_new = new vector<double>(p1.size());
     for ( int j = 0; j < n; ++j )
         (*p_new)[j] = c1 * p1[j] + c2 * p2[j];
     return p_new;
}

void NMSimplex::replace_vertex( int i, vector<double>* p, double fp )
{
    delete vertex_list[i];
    vertex_list[i] = p;
    f_list[i] = fp;
    return;
}

//--------------------------------------------------------------------
