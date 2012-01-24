/** \file gpath_utils.hh
 *  \ingroup libgeom2
 *  \brief Declarations for the C++ geometric-path utitlities.
 *  \author RJG
 *  \version 07-Feb-2008 -- intial code
 *
 */

#ifndef GPATH_UTILS_HH
#define GPATH_UTILS_HH

#include <vector>

#include "../../nm/source/fobject.hh"
#include "geom.hh"
#include "gpath.hh"

XSpline create_quintic_spline(std::vector<Vector3> &P,
			      std::vector<double> &Pdash,
			      double d2a=0.0, double d2b=0.0);
XSpline create_cubic_spline(std::vector<Vector3> &P,
			    std::vector<double> &Pdash,
			    double d2a=0.0, double d2b=0.0);
XPoly create_line(Vector3 &P, double m, double x0, double x1);
int fac(int n);
double choose_func(int n, int i);
double Bernstein(int n, int i, double t);
XBezier create_best_fit_XBezier(std::vector<Vector3> &P, int n, double lslope, double rslope,
			       double lcurv=0.0, double rcurv=0.0);
Bezier create_optimised_Bezier(std::vector<Vector3> &P, int n, double dydxA, double dydxB, double z_pos=0.0);
Bezier create_optimised_Bezier_YZ(std::vector<Vector3> &P, int n, double dydzA, double dydzB, double x_pos=0.0);
Bezier create_optimised_Bezier3D(std::vector<Vector3> &P, int n, double dydxA, double dydxB,
				double dzdxA, double dzdxB);
XPoly create_best_fit_poly(std::vector<Vector3> &P, int n);

class Best_fit_Bezier : public MultivariateFunction {
public:
    Best_fit_Bezier(std::vector<Vector3> points, int order,
		    double dydxA, double dydxB, double z_pos=0.0);
    ~Best_fit_Bezier();
    double eval(std::vector<double> &x);
private:
    std::vector<Vector3> points_;
    int order_;
    double dydxA_;
    double dydxB_;
    double z_pos_;
};

class Best_fit_Bezier_YZ : public MultivariateFunction {
public:
    Best_fit_Bezier_YZ(std::vector<Vector3> points, int order,
		       double dydzA, double dydzB, double x_pos=0.0);
    ~Best_fit_Bezier_YZ();
    double eval(std::vector<double> &x);
private:
    std::vector<Vector3> points_;
    int order_;
    double dydzA_;
    double dydzB_;
    double x_pos_;
};


class Best_fit_Bezier3D : public MultivariateFunction {
public:
    Best_fit_Bezier3D(std::vector<Vector3> points, int order,
		      double dydxA, double dydxB,
		      double dzdxA, double dzdxB);
    ~Best_fit_Bezier3D();
    double eval(std::vector<double> &x);
private:
    std::vector<Vector3> points_;
    int order_;
    double dydxA_;
    double dydxB_;
    double dzdxA_;
    double dzdxB_;
};

int intersect_3D_lines(const Vector3 &P0, const Vector3 &T0,
		       const Vector3 &P2, const Vector3 &T2,
		       double &alf1, double &alf2,
		       Vector3 &P1);

#endif

