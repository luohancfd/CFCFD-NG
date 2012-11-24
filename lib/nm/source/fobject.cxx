/** \file fobject.cxx
 *  \ingroup nm
 *  \brief Implementation of the function-object classes. 
 *  \author PJ
 *  \version 11-Jan-2006 initial coding
 *  \version 19-Oct-2008 Added BivariateFunction (and LinearFunction2)
 *
 * For use in the libgeom2 classes so that we can pass around
 * simple function objects rather than functions themselves.
 *
 */
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "fobject.hh"
#include "../../nm/source/zero_system.hh"
#include "../../nm/source/zero_finders.hh"
#include "../../util/source/useful.h"
using namespace std;

// Base class for functions of one variable.
UnivariateFunction::UnivariateFunction() {}
UnivariateFunction::UnivariateFunction( const UnivariateFunction &f ) {}
UnivariateFunction::~UnivariateFunction() {}
UnivariateFunction* UnivariateFunction::clone() const
{
    return new UnivariateFunction(*this);
}
double UnivariateFunction::eval( double t )
{
    cout << "UnivariateFunction::eval() -- should not be called." << endl;
    return 0.0;
}
string UnivariateFunction::str() const
{
    return "Unknown-Function";
}
void UnivariateFunction::reverse_clustering()
{
    cout << "UnivariateFunction::reverse_clustering() -- should not be called." << endl;
    return;
}

vector<double>* UnivariateFunction::distribute_parameter_values( int n, double t0, double t1 )
{
    vector<double>* tv = new vector<double>(n);
    double dt = 1.0 / (n-1);
    for ( int i = 0; i < n; ++i ) {
	// first mapping is done by a derived function object
	double t = eval(dt * i); 
	(*tv)[i] = (1.0 - t) * t0 + t * t1;  // map to specified range
    }
    return tv;
}

// Overload stream output for UnivariateFunction objects
ostream& operator<<( ostream &os, const UnivariateFunction &f )
{
    os << f.str();
    return os;
}


// Base class for functions of two variables.
BivariateFunction::BivariateFunction() {}
BivariateFunction::BivariateFunction( const BivariateFunction &f ) {}
BivariateFunction::~BivariateFunction() {}
BivariateFunction* BivariateFunction::clone() const
{
    return new BivariateFunction(*this);
}
double BivariateFunction::eval( double r, double s )
{
    cout << "BivariateFunction::eval() -- should not be called." << endl;
    return 0.0;
}
string BivariateFunction::str() const
{
    return "Unknown-BivariateFunction";
}

// Overload stream output for BivariateFunction objects
ostream& operator<<( ostream &os, const BivariateFunction &f )
{
    os << f.str();
    return os;
}


// Base class for functions of several variables.
MultivariateFunction::MultivariateFunction() {}
MultivariateFunction::MultivariateFunction( const MultivariateFunction &f ) {}
MultivariateFunction::~MultivariateFunction() {}
MultivariateFunction* MultivariateFunction::clone() const
{
    return new MultivariateFunction(*this);
}
double MultivariateFunction::eval( vector<double> &x )
{
    cout << "MultivariateFunction::eval() -- should not be called." << endl;
    return 0.0;
}
string MultivariateFunction::str() const
{
    return "Unknown-Function";
}

//--------------------------------------------------------------------------------


// f(t) = m * t + c
LinearFunction::LinearFunction( double _m, double _c )
    : UnivariateFunction(), m(_m), c(_c) {}
LinearFunction::LinearFunction( const LinearFunction &f )
    : UnivariateFunction(), m(f.m), c(f.c) {}
LinearFunction::~LinearFunction() {}
LinearFunction* LinearFunction::clone() const
{
    return new LinearFunction(*this);
}
double LinearFunction::eval( double t )
{
    return m * t + c;
}
string LinearFunction::str() const
{
    ostringstream ost;
    ost << "f(t)=(" << m << " * t + " << c << ")";
    return ost.str();
}
void LinearFunction::reverse_clustering()
{
    // do nothing
}
//--------------------------------------------------------------------------------


// f(t) = y0 * (1-t) + y1 * t
// Actually, we're just going to delegate it all to the original y=mx+c class.
LinearFunction2::LinearFunction2( double y0, double y1 )
    : LinearFunction(y1-y0, y0) {}
LinearFunction2::LinearFunction2( const LinearFunction2 &f )
    : LinearFunction(f.m, f.c) {}
LinearFunction2::~LinearFunction2() {}
//--------------------------------------------------------------------------------


// Roberts' clustering function
RobertsClusterFunction::RobertsClusterFunction( int _end0, int _end1, double _beta )
    : UnivariateFunction(), end0(_end0), end1(_end1), beta(_beta) 
{
    set_underlying_parameters();
}
RobertsClusterFunction::RobertsClusterFunction( const RobertsClusterFunction &f )
    : UnivariateFunction(), end0(f.end0), end1(f.end1), beta(f.beta), 
      alpha(f.alpha), reverse(f.reverse), cluster(f.cluster) {}
RobertsClusterFunction::~RobertsClusterFunction() {}
RobertsClusterFunction* RobertsClusterFunction::clone() const
{
    return new RobertsClusterFunction(*this);
}
double RobertsClusterFunction::eval( double t )
{
    double tbar;
    if ( reverse ) t = 1.0 - t;
    if ( cluster ) tbar = roberts(t); else tbar = t;
    if ( reverse ) tbar = 1.0 - tbar;
    return tbar;
}
string RobertsClusterFunction::str() const
{
    ostringstream ost;
    ost << "RobertsCluster(t)=(end0=" 
	<< end0 << ", end1=" << end1 
	<< ", beta=" << beta << ")";
    return ost.str();
}
void RobertsClusterFunction::reverse_clustering()
{
    int tmp = end0; end0 = end1; end1 = tmp;
    set_underlying_parameters();
    return;
}
void RobertsClusterFunction::set_underlying_parameters()
{
    // Decide on stretching parameters for Robert's transform.
    alpha = 0.5;
    reverse = 0;
    cluster = 1;
    if ( end0 == 0 && end1 == 0 ) cluster = 0;
    if ( beta <= 1.0 ) cluster = 0;
    if ( beta > 1.0 && beta < 1.0001 ) {
	printf( "RobertsClusterFunction: unreasonable beta=%f\n", beta );
	cluster = 0;
    }
    if ( end0 == 1 && end1 == 1 ) alpha = 0.5;
    if ( end0 == 1 && end1 == 0 ) {
	reverse = 1;
	alpha   = 0.0;
    }
    if ( end0 == 0 && end1 == 1 ) {
	reverse = 0;
	alpha   = 0.0;
    }
    return;
}
double RobertsClusterFunction::roberts(double t) const
{
    double lambda, tbar;
    lambda = (beta + 1.0) / (beta - 1.0);
    lambda = pow( lambda, ((t - alpha)/(1.0 - alpha)) );
    tbar = (beta + 2.0 * alpha) * lambda - beta + 2.0 * alpha;
    tbar = tbar / ((2.0 * alpha + 1.0) * (1.0 + lambda));
    return tbar;
}
//--------------------------------------------------------------------------------
/*
Adriaan's implementation of the Valliammai clustering function.
Valliammai, V., Gogoi, A., Grid Quality Improvement using Multi-Block 
Smoothing, Computational Fluid Dynamics Journal, v. 10 no. 2, pp. 169-174, 
July 2001

Author: Adriaan Window
Date: 31-Jan-2006
*/

// constructor
ValliammaiFunction::ValliammaiFunction( double dL0, double dL1, double L, 
					int n )
    : UnivariateFunction(), dL0(dL0), dL1(dL1), L(L), n(n)
{
    cluster = make_cluster(L, n, dL0, dL1, 10, 0.2);
}

ValliammaiFunction::ValliammaiFunction( const ValliammaiFunction &f )
    : UnivariateFunction(), dL0(f.dL0), dL1(f.dL1), L(f.L), n(f.n)
{
    cluster = new vector<double>(f.n);
    for (size_t i=0; i < (*cluster).size(); i++) {
	(*cluster)[i] = (*f.cluster)[i];
    }
}

// destructor
ValliammaiFunction::~ValliammaiFunction() {
    delete cluster;
}

// clone function
ValliammaiFunction* ValliammaiFunction::clone() const
{
    return new ValliammaiFunction(*this);
}

/*
  Function: ValliammaiFunction::eval
  Purpose: Extracts a node point for a given position on the interval.
*/
double ValliammaiFunction::eval( double t )
{
    double alpha, tbar;
    double x, x0, x1, y0, y1;
    int n;
    n = (*cluster).size();
    x = n * t;
    x0 = floor(x);
    x1 = ceil(x);
    if (int(x) != n) {
	y0 = (*cluster)[int(x0)];
	y1 = (*cluster)[int(x1)];
    } else {
	y0 = (*cluster)[n-1];
	y1 = y0;
 	x1 = x0;
 	return 1.0;
    }
    if (x1 == x0){
	alpha = 1.0;
    } else {
	alpha = (x - x0) / (x1 - x0);
    }
    tbar = (1 - alpha) * y0 + alpha * y1;
    return tbar;
}

/*
  Function: ValliammaiFunction::str
  Purpose: Outputs the cluster values to the output stream.
*/
string ValliammaiFunction::str() const
{
    ostringstream ost;
    for (size_t i=0; i < (*cluster).size(); i++) {
	ost << (*cluster)[i] << endl;
    }
    return ost.str();
}

// reverse clustering
void ValliammaiFunction::reverse_clustering()
{
    return;
}


/*
  Function: VallimmaiFunction::valliammai
  Purpose: Evaluate the function denoted as the Valliammai function as per
  reference cited at the beginning of source code.
 */
double ValliammaiFunction::valliammai(double alpha, double t)
{
    return (exp(alpha * t) - 1.0) / (exp(alpha) - 1.0);
}


/*
  Function: VallimmaiFunction::valliammai
  Purpose: Evaluate the function denoted as the Valliammai function as per
  reference cited at the beginning of source code.
 */
vector<double>* ValliammaiFunction::valliammai(double alpha, vector<double>* t)
{
    vector<double>* eta = new vector<double>((*t).size());
    for (size_t i=0; i < (*eta).size(); i++) {
	(*eta)[i] = (exp(alpha * (*t)[i]) - 1.0) / (exp(alpha) - 1.0);
    }
    return eta;
}


/*
  Function: ValliammaiFunction::make_cluster
  Purpose: Formulate the nodal distribution along an interval
*/
vector<double>* ValliammaiFunction::make_cluster( double L, int n,
						  double dL0, double dL1, 
						  int iters, double midFrac)
{
    vector<double> *zeta0, *zeta1, *eta0, *eta1, *eta_final;
    vector<double> *right_piece, *left_piece;
    double alpha, eta_norm;
    int n0, n1;

    /*
      Setup base parameters
    */
    dL0 = dL0 / L * 2.0;
    dL1 = dL1 / L * 2.0;
    n0 = int(0.5 * n);
    n1 = n - n0;


    /*
      Evaluate curve parameter alpha and determine clustering for first
      half of interval.
    */
    zeta0 = new vector<double>(n0);
    for (int i = 0; i < n0; i++) {
	(*zeta0)[i] = i / (n0 - 1.0);
    }
    alpha = find_alpha(dL0, n0);
    eta0 = valliammai(alpha, zeta0);
    for (int i=0; i < n0; i++) {
	(*eta0)[i] = (*eta0)[i] / (*eta0)[n0-1];
    }
    // assign left half
    left_piece = eta0;


    /*
      Evaluate curve parameter alpha and determine clustering for second
      half of interval.
    */
    zeta1 = new vector<double>(n1);
    for (int i = 0; i < n1; i++) {
	(*zeta1)[i] = i / (n1 - 1.0);
    }

    alpha = find_alpha(dL1, n1);
    eta1 = valliammai(alpha, zeta1);
    for (int i = 0; i < n1; i++) {
	(*eta1)[i] = (((*eta1)[i] / (*eta1)[n1-1]) * -1.0) + 2.0;
    }
    // assign right half and flip vector
    right_piece = eta1;
    flip_vector(*right_piece);

    /*
      Concatenate both halfs of interval into single vector
    */
    eta_final = new vector<double>(n-1);
    for (size_t i=0; i < (*left_piece).size(); i++) {
      (*eta_final)[i] = (*left_piece)[i];// / *((*right_piece).end()-1);
    }
    for (size_t i=(*left_piece).size(); i < ((*left_piece).size()+
					  (*right_piece).size()-1); i++) {
	(*eta_final)[i] = (*right_piece)[i-(*left_piece).size()];
    }

    /*
      Check end values and enforce 0.0 to 1.0 policy.
    */
    if ((*eta_final)[0] != 0.0) {
	(*eta_final)[0] = 0.0;
    }

    /*
      Normalise all values to 1.0
    */
    eta_norm = *((*eta_final).end()-1);
    cout << "norm = " << eta_norm << endl;
    for (size_t i=0; i < (*eta_final).size(); i++) {
        (*eta_final)[i] = (*eta_final)[i]/eta_norm;
    }
    // ditch last value
    (*eta_final).pop_back();

    /*
      Perform curve smoothing operation on mid portion of cluster.
      eta_final is modified within the function.
    */
    smooth_curve(*eta_final, midFrac, iters);
    delete zeta0;
    delete zeta1;
    delete eta0;
    delete eta1;
    
    return eta_final;
}


/*
  Function: ValliammaiFunction::flip_vector
  Purpose: Flip the elements in a vector<>
*/
void ValliammaiFunction::flip_vector(vector<double>& v)
{
    vector<double>* v_star = new vector<double>(v.size());
    for (size_t i = 0; i < v.size(); i++) {
	(*v_star)[i] = v[(v.size()-1)-i];
    }
    v = *v_star;
    delete v_star;
}


/*
  Function: ValliammaiFunction::find_alpha
  Purpose: Determine the value of alpha require for the desired clustering.
*/
double ValliammaiFunction::find_alpha(double dL, int n)
{
    /*
      Define logical space finite length scale and initial vector
    */
    double dx, da;
    double nominal_alpha;
    vector<double> *x, *alphav, *s, *ssorted;
    x = new vector<double>(n);
    nominal_alpha = 1.0;
    dx = 1.0 / (n - 1);
    for (size_t i=0; i < (*x).size(); i++) {
	(*x)[i] = i * dx;
    }
    
    int alpha_steps = int(floor((10.0 - 1.0) / 0.0001));
    da = (10.0 - 1.0) / alpha_steps;
    alphav = new vector<double>(alpha_steps);
    s = new vector<double>(alpha_steps);
    for (int i=0; i < alpha_steps; i++) {
	(*alphav)[i] = 1.0 + i * da;
    }
    for (int i=0; i < alpha_steps; i++) {
	(*s)[i] = valliammai( (*alphav)[i], (*x)[1] ) 
	    - valliammai( (*alphav)[i], (*x)[0] );
    }

    // find maximum spacing
    ssorted = new vector<double>((*s).size());
    copy((*s).begin(), (*s).end(), (*ssorted).begin());
    sort( (*ssorted).begin(), (*ssorted).end() );
    vector<double>::iterator s_max = max_element( (*s).begin(), (*s).end() );
    vector<double>::iterator space = find( (*s).begin(), (*s).end(), 
	*(lower_bound( (*ssorted).begin(), (*ssorted).end(), dL )) );
    
    int s_idx;
    if (*s_max < dL) {
	s_idx = distance( (*s).begin(), s_max );
	dL = *s_max;
	nominal_alpha = (*alphav)[s_idx];
    }
    else {
	s_idx = distance( (*s).begin(), space );
	nominal_alpha = (*alphav)[s_idx];
    }
    delete alphav;
    delete ssorted;
    delete s;
    delete x;
    return nominal_alpha;
}


/*
  Function: ValliammaiFunction::smooth_curve
  Purpose: This function performs a smoothing operation on a portion of
  the curve passed in. The vector is modified within this function and
  so nothing is returned.
 */
void ValliammaiFunction::smooth_curve(vector<double>& t, double midFrac, 
				      int iters)
{
    /*
      Define and extract portion of curve to be smoothed.
     */
    vector<double> eta_mid;
    int mid_start, mid_end;
    mid_start = int(floor(t.size() / 2 
		      - midFrac / 2 * t.size()));
    mid_end = int(floor(t.size()/2 
		    + midFrac / 2 * t.size()));
    for (int i=mid_start; i < mid_end; i++) {
	(eta_mid).push_back(t[i]);
    }
    /*
      Perform smoothing
    */
    vector<double> eta_mid_star = eta_mid;
    for (int count=0; count < iters; count++) {
	for (size_t i=1; i < (eta_mid).size()-1; i++) {
	    (eta_mid_star)[i] = ( (eta_mid)[i+1] + (eta_mid)[i-1]) / 2.0;
	}
	eta_mid = eta_mid_star;
    }
    /*
      Replace smoothed values in original vector.
    */
    for (int i=mid_start; i < mid_end; i++) {
	t[i] = (eta_mid)[i-mid_start];
    }
}

//--------------------------------------------------------------------------------


// Combination of two UnivariateFunction's
// gamma: location in t space of the discontinuity
// uf0: univariate function used for 0 < t < gamma
// uf1: univariate function used for gamma < t < 1
DiscontinuousUnivariateFunction::
DiscontinuousUnivariateFunction( double _gamma, UnivariateFunction * _uf0, UnivariateFunction * _uf1 )
    : UnivariateFunction(), gamma( _gamma )
{
    uf0 = _uf0->clone();
    uf1 = _uf1->clone();
}
DiscontinuousUnivariateFunction::
DiscontinuousUnivariateFunction( const DiscontinuousUnivariateFunction &f )
    : UnivariateFunction(), gamma( f.gamma )
{
    uf0 = f.uf0->clone();
    uf1 = f.uf1->clone();
}
DiscontinuousUnivariateFunction::
~DiscontinuousUnivariateFunction()
{
    delete uf0;
    delete uf1;
}
DiscontinuousUnivariateFunction*
DiscontinuousUnivariateFunction::clone() const
{
    return new DiscontinuousUnivariateFunction(*this);
}
double DiscontinuousUnivariateFunction::eval( double t )
{
    double tbar;

    if ( t < gamma )
   	tbar = uf0->eval( t/gamma ) * gamma;
    else
        tbar = uf1->eval( ( t - gamma ) / ( 1.0 - gamma ) ) * ( 1.0 - gamma ) + gamma;
    
    return tbar;
}
string DiscontinuousUnivariateFunction::str() const
{
    ostringstream ost;
    ost << "Function 0: " << uf0->str() << endl;
    ost << "Function 1: " << uf1->str() << endl;
    return ost.str();
}
void DiscontinuousUnivariateFunction::reverse_clustering()
{
    gamma = 1.0 - gamma;
    uf0->reverse_clustering();
    uf1->reverse_clustering();
    
    return;
}

//--------------------------------------------------------------------



// f(r,s) = area-weighting of corner values
BilinearFunction::BilinearFunction( double _v0, double _v1, double _v2, double _v3 )
    : BivariateFunction(), v0(_v0), v1(_v1), v2(_v2), v3(_v3) {}
BilinearFunction::BilinearFunction( const BilinearFunction &f )
    : BivariateFunction(), v0(f.v0), v1(f.v1), v2(f.v2), v3(f.v3) {}
BilinearFunction::~BilinearFunction() {}
BilinearFunction* BilinearFunction::clone() const
{
    return new BilinearFunction(*this);
}
double BilinearFunction::eval( double r, double s )
{
    return v0*(1.0-r)*(1.0-s) + v1*r*(1.0-s) + v2*r*s + v3*(1.0-r)*s ;
}
string BilinearFunction::str() const
{
    ostringstream ost;
    ost << "f(r,s)=Area-weighting-of(v0=" << v0 << ", v1=" << v1 
	<< ", v2=" << v2 << ", v3=" << v3 << ")";
    return ost.str();
}

//--------------------------------------------------------------------
// Simple Hyperbolic-tangent clustering
// Ref: http://www.cfd-online.com/Wiki/Structured_mesh_generation

class VinokurFunction : public ZeroFunction  {
public:
    VinokurFunction( double _B );
    VinokurFunction(const VinokurFunction& vfun);
    virtual ~VinokurFunction();
    virtual int eval(double x, double &y); // y = f(x)
private:
    double B;
};

VinokurFunction::VinokurFunction( double _B ) : B( _B ) {}
VinokurFunction::VinokurFunction( const VinokurFunction &z )
: B( z.B ) {}
VinokurFunction::~VinokurFunction() {}
int VinokurFunction::eval( double x, double &y)
{
    y = B - sinh(x) / x;

    return SUCCESS;
}

// Hyberbolic-tangent clustering function
HypertanClusterFunction::HypertanClusterFunction( double _dL0, double _dL1 )
    : UnivariateFunction(), dL0(_dL0), dL1(_dL1)
{
    set_underlying_parameters();
}
HypertanClusterFunction::HypertanClusterFunction( const HypertanClusterFunction &f )
    : UnivariateFunction(), dL0(f.dL0), dL1(f.dL1)
{
    set_underlying_parameters();
}
HypertanClusterFunction::~HypertanClusterFunction() {}
HypertanClusterFunction* HypertanClusterFunction::clone() const
{
    return new HypertanClusterFunction(*this);
}
double HypertanClusterFunction::eval( double t )
{
    double u = 0.5 * ( 1.0 + tanh( delta * ( t - 0.5 ) ) );
    double tbar = u / ( A + ( 1.0 - A ) * u );
    return (tbar-tbar0)*tbar_scale;
}
string HypertanClusterFunction::str() const
{
    ostringstream ost;
    ost << "HypertanClusterFunction(t)=(dL0="
        << dL0 << ", dL1=" << dL1 << ")";
    return ost.str();
}
void HypertanClusterFunction::reverse_clustering()
{
    int tmp = dL0; dL0 = dL1; dL1 = tmp;
    set_underlying_parameters();
    return;
}
void HypertanClusterFunction::set_underlying_parameters()
{
    // Solve for delta
    A = sqrt( dL1 ) / sqrt( dL0 );
    B = 1.0 / sqrt( dL1 * dL0 );
    VinokurFunction f = VinokurFunction( B );
    Bisection bsm = Bisection( &f, 1.0e-12 );
    delta = bsm( 1.0e-10, 1.0e3 );
    // Solve for shifting parameters
    tbar0 = 0.0; tbar_scale = 1.0;
    tbar0 = this->eval(0);
    tbar_scale = 1.0 / this->eval(1.0);
    // cout << "tbar0 = " << tbar0 << ", tbar_scale = " << tbar_scale << endl;

    return;
}
