
#ifndef BEST_FIT_SPLINE_HH
#define BEST_FIT_SPLINE_HH

#include <string>
#include <vector>
#include <valarray>

#include "../../nm/source/fobject.hh"
#include "geom.hh"


class Best_fit_spline : public MultivariateFunction {
public:
    Best_fit_spline(std::string input_file,
		    double errwt1, double errwt2,
		    double y0, double x1);
    Best_fit_spline(const Best_fit_spline &bs);
    ~Best_fit_spline();
    Best_fit_spline* clone() const;
    double eval( std::vector<double> &x );
    std::string str() const;
private:
    double errwt1_;
    double errwt2_;
    double y0_;
    double x1_;
    std::vector<Vector3> points_;
};

#endif
