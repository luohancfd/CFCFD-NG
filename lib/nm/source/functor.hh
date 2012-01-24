// Author: Rowan J. Gollan
// Date: 08-Jul-2008

#ifndef FUNCTOR_HH
#define FUNCTOR_HH

class Univariate_functor {
public:
    virtual double operator()(double x) = 0;
    virtual ~Univariate_functor() {}
    virtual Univariate_functor* clone() const
    { return 0; }
};

class Bivariate_functor {
public:
    virtual double operator()(double x, double y) = 0;
    virtual ~Bivariate_functor() {}
    virtual Bivariate_functor* clone() const
    { return 0; }
};

#endif
