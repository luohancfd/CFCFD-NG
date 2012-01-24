/** \file fobject.hh
 *  \ingroup nm
 *  \brief Declarations for function objects for use in libgeom2.
 *  \author PJ
 *  \version 11-Jan-2006
 *
 */

#ifndef FOBJECT_HH
#define FOBJECT_HH

#include <string>
#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

// Base class for functions of one variable.
class UnivariateFunction {
public:
    UnivariateFunction();
    UnivariateFunction( const UnivariateFunction &f );
    virtual ~UnivariateFunction();
    virtual UnivariateFunction* clone() const;
    virtual double eval( double t );
    virtual string str() const;
    // reverse_clustering() alters the internal parameters so that 
    // the clustering (i.e. deviation from a linear distribution)
    // over the t-parameter range 0.0-->1.0 looks like
    // the original evluated over the range 1.0-->0.0
    // We still want the parameter to vary from 0.0-->1.0.
    virtual void reverse_clustering(); 
    // Utility function for grid generation.
    vector<double>* distribute_parameter_values( int n, double t0=0.0, double t1=1.0 );
};

// Overload stream output for UnivariateFunction objects
ostream& operator<<( ostream &os, const UnivariateFunction &f );


// Base class for functions of two variables.
class BivariateFunction {
public:
    BivariateFunction();
    BivariateFunction( const BivariateFunction &f );
    virtual ~BivariateFunction();
    virtual BivariateFunction* clone() const;
    virtual double eval( double r, double s );
    virtual string str() const;
};

// Overload stream output for BivariateFunction objects
ostream& operator<<( ostream &os, const BivariateFunction &f );


// Base class for functions of several variables.
class MultivariateFunction {
public:
    MultivariateFunction();
    MultivariateFunction( const MultivariateFunction &f );
    virtual ~MultivariateFunction();
    virtual MultivariateFunction* clone() const;
    virtual double eval( vector<double> &x );
    virtual string str() const;
};

//-----------------------------------------------------------------------------
// y = m * x + c
class LinearFunction : public UnivariateFunction {
public:
    double m; // slope
    double c; // intercept
    LinearFunction( double _m=1.0, double _c=0.0 );
    LinearFunction( const LinearFunction &f );
    virtual ~LinearFunction();
    virtual LinearFunction* clone() const;
    virtual double eval( double t );
    virtual string str() const;
    virtual void reverse_clustering(); 
};

// y = y0*(1-x) + y1*x
class LinearFunction2 : public LinearFunction {
public:
    LinearFunction2( double y0=0.0, double y1=1.0 );
    LinearFunction2( const LinearFunction2 &f );
    virtual ~LinearFunction2();
};


// Roberts' clustering function
class RobertsClusterFunction : public UnivariateFunction {
public:
    bool   end0, end1; // flags to indicate clustering to each end 
    double beta; // stretching factor 1.0 < beta < +inf
                 // closer to 1.0 gives stronger clustering
    RobertsClusterFunction( int _end0=0, int _end1=0, double _beta=0.0 );
    RobertsClusterFunction( const RobertsClusterFunction &f );
    virtual ~RobertsClusterFunction();
    virtual RobertsClusterFunction* clone() const;
    virtual double eval( double t );
    virtual string str() const;
    virtual void reverse_clustering(); 
private:
    // underlying parameters...
    double alpha; // location of stretching: 
                  // 0.5 cluster at both ends
                  // 0.0 cluster toward t=1.0
    bool   reverse;
    bool   cluster;
    void set_underlying_parameters();
    double roberts(double t) const; // returns transformed ordinate
};


// valliammai clustering
class ValliammaiFunction : public UnivariateFunction {
public:
    double dL0, dL1; // size of cell at each end
    double L; // overall length of the path to be clustered
    int n; // number of nodes to be distributed along path
    
    ValliammaiFunction( double dL0, double dL1, double L, int n);
    ValliammaiFunction( const ValliammaiFunction &f );
    virtual ~ValliammaiFunction();
    virtual ValliammaiFunction* clone() const;
    virtual double eval( double t );
    virtual string str() const;
    virtual void reverse_clustering();
    vector<double>* cluster;
private:
    vector<double>* make_cluster(double L, int n, double dL0, double dL1,
				 int iters=10, double midFrac=0.2);
    double find_alpha(double dL, int n);
    void smooth_curve(vector<double>&, double, int);
    double valliammai(double alpha, double t);
    vector<double>* valliammai(double alpha, vector<double>*);
    void flip_vector(vector<double>&);
};    

// Combination of two UnivariateFunction's
class DiscontinuousUnivariateFunction : public UnivariateFunction {
public:
    bool   end0, end1; // flags to indicate clustering to each end 
    double beta; // stretching factor 1.0 < beta < +inf
                 // closer to 1.0 gives stronger clustering
    DiscontinuousUnivariateFunction( double gamma, UnivariateFunction * _uf0, UnivariateFunction * _uf1 );
    DiscontinuousUnivariateFunction( const DiscontinuousUnivariateFunction &f );
    virtual ~DiscontinuousUnivariateFunction();
    virtual DiscontinuousUnivariateFunction* clone() const;
    virtual double eval( double t );
    virtual string str() const;
    virtual void reverse_clustering(); 
private:
    UnivariateFunction * uf0;
    UnivariateFunction * uf1;
    double gamma;
};

// Adriaan's Hyperbolic-tangent clustering...
// class HyptanClusterFunction : public UnivariateFunction {
// public:
//     double dL0, dL1; // size of cell at each end  
//     double L; // overall length of the Path along which the points are distributed
//     int n; // number of points to be spread along the path (including end points)
//     HyptanClusterFunction( double dL0, double dL1, double L=1.0, int n=10 );
//     HyptanClusterFunction( const HyptanClusterFunction &f );
//     virtual ~HyptanClusterFunction();
//     virtual HyptanClusterFunction* clone() const;
//     virtual double eval( double t ) const;
//     virtual string str() const;
//     virtual void reverse_clustering(); 
// private:
//     double beta; // any bits that are needed during the calculation...
// };


//-----------------------------------------------------------------------------
// f(r,s) = Area-weighting of corner values.
//
//     1    3-------2
//     |    |       |
//     s    |       |
//     |    |       |
//     0    0-------1
//
//          0-- r --1
class BilinearFunction : public BivariateFunction {
public:
    double v0, v1, v2, v3;
    BilinearFunction( double _v0, double _v1, double _v2, double _v3 );
    BilinearFunction( const BilinearFunction &f );
    virtual ~BilinearFunction();
    virtual BilinearFunction* clone() const;
    virtual double eval( double r, double s );
    virtual string str() const;
};

#endif
