// Author: Rowan J. Gollan
// Date: 10-Sep-2008

#ifndef NURBS_HH
#define NURBS_HH

#include <vector>
#include "geom.hh"

int find_span(double u, int n, int p, const std::vector<double> &U);
int basis_funs(double u, int i, int p,
	       const std::vector<double> &U, std::vector<double> &N);
double one_basis_fun(int p, const std::vector<double> U,
		     int i, double u);
class Mapped_point {
public:
    Mapped_point();
    Mapped_point(double wx, double wy, double wz, double w);
    Mapped_point(const Mapped_point &mp);

    ~Mapped_point();

    Mapped_point& operator=(const Mapped_point &mp);

    Vector3 to_Vector3();

    // I normally wouldn't expose these members as public
    // but this is a very simple class.  Orginally it was
    // a data structure, but it became a class so that I
    // could define the multiply by scalar operator.

    double wx;
    double wy;
    double wz;
    double w;
};

Mapped_point operator+(const Mapped_point &m1, const Mapped_point &m2);
Mapped_point operator-(const Mapped_point &m1, const Mapped_point &m2);
Mapped_point operator*(double scalar, const Mapped_point &mp);
Mapped_point operator*(const Mapped_point &mp, double scalar);

double abs(const Mapped_point &m);

Vector3 nurbs_curve_point(double u, int p,
			  const std::vector<double> &U,
			  const std::vector<Mapped_point> &Pw);

Vector3 nurbs_surface_point(double u, int p,
			    const std::vector<double> &U,
			    double v, int q,
			    const std::vector<double> &V,
			    const std::vector<std::vector<Mapped_point> > &Pw);

#endif
