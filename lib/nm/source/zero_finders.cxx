/** \file zero_finders.cxx
 * \ingroup nm
 * \brief Class and method definitions for the collection of
 *        zero-finders.
 *
 * 
 * \author Rowan J Gollan
 * \version 21-Apr-2006 - initial coding
 * \version 26-Nov-2013 - change to using vector instead of valarray (PJ)
 *
 **/

#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../util/source/useful.h"
#include "no_fuss_linear_algebra.hh"
#include "zero_system.hh"
#include "zero_finders.hh"

using namespace std;

/// \brief Default constructor
///
/// \author Rowan J. Gollan
/// \version 07-Jan-07
///
ZeroFinder::ZeroFinder() {}

/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 21-Apr-2006
///
ZeroFinder::ZeroFinder( int ndim, double tolerance, int max_iterations )
    : ndim_( ndim ), tol_( tolerance ), max_iter_( max_iterations ) {}

/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 21-Apr-2006
///
ZeroFinder::ZeroFinder( const ZeroFinder &z )
    : ndim_( z.ndim_ ), tol_( z.tol_ ), max_iter_( z.max_iter_ ) {}

/// \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 21-Apr-2006
///
ZeroFinder::~ZeroFinder() {}

/// \brief Set constants for the zero finder.
///
/// \author Rowan J. Gollan
/// \version 07-Jan-07
///
void ZeroFinder::set_constants(int ndim, double tolerance, int max_iterations)
{
    ndim_ = ndim;
    tol_ = tolerance;
    max_iter_ = max_iterations;
}

/// \brief  Bisection method functor
/// \author Brendan T. O'Flaherty

Bisection::
Bisection() {}

Bisection::
Bisection(const Bisection &b)
{
    tol_ = b.tol_;
    zfun_ = b.zfun_;
}

// Bisection::
// Bisection(double tol)
// {
//     tol_ = tol;
// }

Bisection::
Bisection(ZeroFunction *zfun,
	  double tol)
{
    tol_ = tol;
    zfun_ = zfun;
}

Bisection::
~Bisection() {}

Bisection*
Bisection::
clone() const {
    return new Bisection(*this);
}

double
Bisection::
operator()(double x_left, 
	   double x_right)
{
    if (x_left > x_right) {
	printf("error in bisection root:\n");
	printf("x_lower > x_upper\n");
	printf("swap bounds and try again.\n");
	exit(1);
    }
    double y_left;
    double y_right;
    if (zfun_->eval(x_left, y_left) != SUCCESS) { 
	printf("error in bisection root:\n");
	printf(" failure to evaluate left bound, x = %g\n", x_left);
	printf(" exiting to system...\n");
	exit(1);
    }
    if (zfun_->eval(x_right, y_right) != SUCCESS) {
	printf("error in bisection root:\n");
	printf(" failure to evaluate right bound, x = %g\n", x_right);
	printf(" exiting to system...\n");
	exit(1);
    }
    if (y_left*y_right >= 0.0) {
	printf("error in bisection root:\n");
	printf("x-bounds [%g, %g] correspond to y-bounds [%g, %g]\n", 
	       x_left, x_right, y_left, y_right);
	printf("which do not fall either side of zero.\n");
	exit(1);
    }
    
    double x_mid;
    // We refine the x-bracket until it is small enough.
    for(x_mid=(x_left+x_right)/2.0; fabs(x_left-x_mid) > tol_; x_mid=(x_left+x_right)/2.0) {
	if (zfun_->eval(x_left, y_left) != SUCCESS) { 
	    printf("error in bisection root:\n");
	    printf(" failure to evaluate left bound, x = %g\n", x_left);
	    printf(" exiting to system...\n");
	    exit(1);
	}
	if (zfun_->eval(x_mid, y_right) != SUCCESS) {
	    printf("error in bisection root:\n");
	    printf(" failure to evaluate right bound, x = %g\n", x_mid);
	    printf(" exiting to system...\n");
	    exit(1);
	}
	if (y_left*y_right <= 0.0) {
	    x_right = x_mid; // use left interval
	} else {
	    x_left = x_mid; // use right interval
	}
    }
    return x_mid;
}

/// \brief  Muller's method
/// \author Brendan T. O'Flaherty

Muller::
Muller() {}

Muller::
Muller(const Muller &b)
{
    tol_ = b.tol_;
    zfun_ = b.zfun_;
}

// Muller::
// Muller(double tol)
// {
//     tol_ = tol;
// }

Muller::
Muller(ZeroFunction *zfun,
       double tol)
{
    tol_ = tol;
    zfun_ = zfun;
}

Muller::
~Muller() {}

Muller*
Muller::
clone() const {
    return new Muller(*this);
}

double
Muller::
operator()(double x_left, 
	   double x_right)
{
    if (x_left > x_right) {
	printf("error in muller root:\n");
	printf("x_lower > x_upper\n");
	printf("%5.4e > %5.4e\n", x_left, x_right);
	printf("swap bounds and try again.\n");
	exit(1);
    }
    int iter = 0;
    int ITER_MAX = 100;
    double y_left;
    double y_right;
    if (zfun_->eval(x_left, y_left) != SUCCESS) { 
	printf("error in muller root:\n");
	printf(" failure to evaluate left bound, x = %g\n", x_left);
	printf(" exiting to system...\n");
	exit(1);
    }
    if (zfun_->eval(x_right, y_right) != SUCCESS) {
	printf("error in muller root:\n");
	printf(" failure to evaluate right bound, x = %g\n", x_right);
	printf(" exiting to system...\n");
	exit(1);
    }
    if (y_left*y_right >= 0.0) {
	printf("error in muller root:\n");
	printf("x-bounds [%g, %g] correspond to y-bounds [%g, %g]\n", 
	       x_left, x_right, y_left, y_right);
	printf("which do not fall either side of zero.\n");
	exit(1);
    }

    // at this point we have checked the bounds. 
    // this requires at minimum two evaluations.

    // Muller's method, p278, Press W.H. et al. (1988). "Numerical Recipes in C". Cambridge University Press.
    // like the bisection method, this method requires two evaluations per iteration
    // however it uses a polynomial curve fit to find a better approximation of the root
    // and so converges faster (for a well behaved curve).
    
    double q,A,B,C,den,x_new1,x_new2;
    double x_root = (x_left+x_right)/2.0; // initial guess
    double y_root = 1; // ensure at least one evaluation
    
    for (x_mid_=x_root; fabs(y_root) > tol_; x_mid_=(x_left+x_right)/2.0) {
	if (zfun_->eval(x_mid_, y_mid_) != SUCCESS) {
	    printf("error in muller root:\n");
	    printf(" failure to evaluate right bound, x = %g\n", x_mid_);
	    printf(" exiting to system...\n");
	    exit(1);
	}
	
        q = (x_right - x_mid_)/(x_mid_ - x_left);
	A = q*y_right - q*(1 + q)*y_mid_ + q*q*y_left;
	B = (2*q + 1)*y_right - (1 + q)*(1 + q)*y_mid_ + q*q*y_left;
	C = (1 + q)*y_right;
	
	den =  B*B - 4*A*C;
	if (den < 0) {
	    printf("error in muller root:\n");
	    printf(" complex root when evaluating polynomial.");
	    printf(" exiting to system...\n");
	    exit(1);
	}
	
	x_new1 = x_right - (x_right - x_mid_)*(2*C/(B + pow(den,0.5)));
	x_new2 = x_right - (x_right - x_mid_)*(2*C/(B - pow(den,0.5)));
	
	// now we have five points
	// the original bounds, the midpoint, and the two polynomial roots
	// first check which root is within the original bounds
	
	if ((x_new1 > x_left) && (x_new1 < x_right)) {
	    x_root = x_new1;
	} else {
	    x_root = x_new2;
	}
	
	// now we have four points
	// the original bounds, the midpoint, and a good approximation.
	// find the new bounds.
	
	zfun_->eval(x_root, y_root);
	
	if (y_left*y_mid_ <= 0.0) { // use left interval
	    x_right = x_mid_;
	    y_right = y_mid_;
	    if (y_root*y_right <= 0.0) { // use right interval
		x_left = x_root;
		y_left = y_root;
	    } else { // use the left interval
		x_right = x_root;
		y_right = y_root;
	    }
	} else { // use right interval
	    x_left = x_mid_;
	    y_left = y_mid_;
	    if (y_root*y_right <= 0.0) { // use right interval
		x_left = x_root;
		y_left = y_root;
	    } else { // use left interval
		x_right = x_root;
		y_right = y_root;
	    }
	}
	++iter;
	if (iter >= ITER_MAX) {
	    printf("Error in Muller zero method\n");
	    printf("Max number of iterations reached (%i)\n", ITER_MAX);
	    printf("Current bounds: [%4.3e, %4.3e], [%4.3e, %4.3e]\n", x_left, x_right, y_left, y_right);
	    printf("Exiting to system...\n");
	}
	//printf("%g %g, %g %g\n", x_left, x_right, y_left, y_right);
    }
    
    return x_root;
}

/// \brief Default constructor
///
/// \author Rowan J. Gollan
/// \version 07-Jan-07
///
NewtonRaphsonZF::NewtonRaphsonZF()
    : ZeroFinder() {}

/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 21-Apr-2006
///
NewtonRaphsonZF::NewtonRaphsonZF( int ndim,
				  double tolerance,
				  int max_iterations,
				  bool has_Jacobian )
    : ZeroFinder( ndim, tolerance, max_iterations ), has_Jac_( has_Jacobian )
{
    G_.resize( ndim );
    minusG_.resize( ndim );
    dGdy_.resize( ndim, ndim );
    y_old_.resize( ndim );
    y_new_.resize( ndim );
    dely_.resize( ndim );

}

/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 21-Apr-2006
///
NewtonRaphsonZF::NewtonRaphsonZF( const NewtonRaphsonZF &n )
    : ZeroFinder( n.ndim_, n.tol_, n.max_iter_ ), has_Jac_( n.has_Jac_ )
{
    G_.resize( ndim_ );
    minusG_.resize( ndim_ );
    dGdy_.resize( ndim_, ndim_ );
    y_old_.resize( ndim_ );
    y_new_.resize( ndim_ );
    dely_.resize( ndim_ );
}

/// \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 21-Apr-2006
///
NewtonRaphsonZF::~NewtonRaphsonZF() {}


/// \brief Set constants for zero finder.
///
/// \author Rowan J. Gollan
/// \version 07-Jan-07
///
void NewtonRaphsonZF::set_constants(int ndim, double tolerance,
				   int max_iterations, bool has_Jacobian)
{
    ZeroFinder::set_constants(ndim, tolerance, max_iterations);
    has_Jac_ = has_Jacobian;
    G_.resize( ndim );
    minusG_.resize( ndim );
    dGdy_.resize( ndim, ndim );
    y_old_.resize( ndim );
    y_new_.resize( ndim );
    dely_.resize( ndim );

}

/// \brief Newton-Raphson solving algorithm
///
/// \author Rowan J Gollan
/// \version 21-Apr-2006
///
int NewtonRaphsonZF::solve( ZeroSystem &zsys, const vector<double> &y_guess, vector<double> &y_out )
{
    int count;

    zsys.f( y_guess, G_ );
    scale_vector2vector( G_, -1.0, minusG_);
    if( has_Jac_ ) {
	zsys.Jac( y_guess, dGdy_ );
    }
    else {
	;// We need some derivative technique
    }
    copy_vector( y_guess, y_old_ );

    for( count = 0; count < max_iter_; ++count ) {
	gaussian_elimination( dGdy_, dely_, minusG_ );
	add_vectors( y_new_, y_old_, dely_ );
	// cout << "count = " << count << ", y_new_ = "; print_vector(y_new_);
	// cout << "dely = "; print_vector(dely_);
	if ( test_tol() )
	    break;
	else {
	    copy_vector( y_new_, y_old_ );
	    zsys.f( y_old_, G_ );
	    scale_vector2vector( G_, -1.0, minusG_);
	    if( has_Jac_ ) {
		zsys.Jac( y_old_, dGdy_ );
	    }
	    else {
		;// We need some derivative technique
	    }
	}
    }

    if ( count >= max_iter_ ) {
	cerr << "NewtonRaphsonZF::solve()\n";
	cerr << "Newton-Raphson method did not coverge (in zero_finders.cxx)\n";
	cerr << "Iterations = " << count << ", tolerance= " << tol_ << endl;
	cerr << "Guessed input= \n";
	print_vector(y_guess);
	cerr << "Values after iteration= \n";
	print_vector(y_new_);
	cerr << "dely values after iteration= \n";
	print_vector(dely_);
	return FAILURE;
    }

    for( size_t i = 0; i < y_out.size(); ++i ) y_out[i] = y_new_[i];

    return SUCCESS;
}


bool NewtonRaphsonZF::test_tol()
{
    bool flag = true;
    for( size_t i = 0; i < dely_.size(); ++i ) {
	if( fabs(dely_[i]) >= tol_ ) {
	    flag = false;
	}
    }

    return flag;
}


