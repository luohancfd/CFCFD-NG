// Author: Rowan J. Gollan
// Date: 04-Nov-2008
// Place: Hampton, Virginia, USA

#ifndef SEGMENTED_FUNCTOR_HH
#define SEGMENTED_FUNCTOR_HH

#include <vector>

#include "functor.hh"

class Segmented_functor : public Univariate_functor {
public:
    Segmented_functor(std::vector<Univariate_functor*> &f, std::vector<double> &breaks);
    Segmented_functor(const Segmented_functor &s);
    ~Segmented_functor();
    Segmented_functor* clone() const;
    
    double operator()(double x);
    double get_break(int i) { return breaks_[i]; }
    double size() { return breaks_.size(); }

private:
    std::vector<Univariate_functor*> f_;
    std::vector<double> breaks_;
};

#endif
