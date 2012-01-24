// Author: Rowan J. Gollan
// Date: 02-Mar-2010
// Place: Poquoson, Virginia, USA
//
// This code is adapted from...
//
// Chapter 2.3 of
// Press, W.H., Teukolsky, S.A., Vetterling, W.T and Flannery, B.P. (2007)
// Numerical Recipes: The Art of Scientific Computing, 3rd Edition,
// Cambridge University Press, New York
//

#include <vector>
#include "no_fuss_linear_algebra.hh"

class LUdcmp {
public:
    LUdcmp(const Valmatrix &a);
    void solve(const std::vector<double> &b, std::vector<double> &x);
    Valmatrix lu()
    { return lu_; }
private:
    size_t n_;
    Valmatrix lu_;
    std::vector<size_t> indx_;
    double d_;
};

