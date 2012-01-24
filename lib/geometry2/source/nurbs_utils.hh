// Author: Rowan J. Gollan
// Date: 19-May-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#ifndef NURBS_UTIL_HH
#define NURBS_UTIL_HH

#include <vector>

#include "nurbs.hh"
#include "geom.hh"
#include "gpath.hh"
#include "surface.hh"

int curve_knot_ins(int np, int p, const std::vector<double> &UP,
		   const std::vector<Mapped_point> &Pw, double u, 
		   int k, int s, int r, int &nq, 
		   std::vector<double> &UQ, std::vector<Mapped_point> &Qw);
int refine_knot_vector_curve(int n, int p, const std::vector<double> &U,
			     const std::vector<Mapped_point> &Pw,
			     const std::vector<double> &X, int r,
			     std::vector<double> &Ubar,
			     std::vector<Mapped_point> &Qw);
int estimate_tangents(const std::vector<Vector3> &Q, std::vector<Vector3> &T,
		      bool preserve_corners);
Nurbs local_rat_quad_curve_interp(const std::vector<Vector3> &Q, bool preserve_corners=true);
int quad_curve_interp(const Vector3 &Q0, const Vector3 &T0,
		      const Vector3 &Q1, const Vector3 &T1,
		      std::vector<Vector3> &R,
		      std::vector<double> &w);


int quad_curve_interp_weight(const Vector3 &Q0, const Vector3 &R1, const Vector3 &Q1,
			     double &w);
int make_one_arc(const Vector3 &P0, const Vector3 &T0, const Vector3 &P2, const Vector3 &T2,
		 const Vector3 &P, Vector3 &P1, double &w1);

void deriv_basis_funs(int i, double u, int p, int n,
		      const std::vector<double> &U, std::vector<std::vector<double> > &ders);

void curve_derivs(int n, int p, const std::vector<double> &U,
		  const std::vector<Mapped_point> &Pw,
		  double u, int d,
		  std::vector<Mapped_point> &CKw);
int factorial(int n);
int bin_coeff(int n, int k);
void rat_curve_derivs(const std::vector<Vector3> &Aders, const std::vector<double> &wders,
		      int d, std::vector<Vector3> &CK);
void rat_curve_derivs(const std::vector<Mapped_point> &CKw, int d, std::vector<Vector3> &CK);
double dist_point_projection(const Vector3 &P, const Nurbs &C,
			     double &t_found, Vector3 &C_found,
			     double eps1=1.0e-6, double eps2=1.0e-8);

int fit_with_conic(int ks, int ke, const std::vector<Vector3> &Q,
		   const Vector3 &Ts, const Vector3 &Te,
		   double E, Mapped_point &Rw);

int search_max_span(int ks, int ke,
		    const std::vector<Vector3> &Q, const std::vector<Vector3> &T,
		    double E);
Nurbs quad_curve_approx(const std::vector<Vector3> &Q, double E, int Kmax=-1,
			bool preserve_corners=true);
int interp_homogeneous_points(const std::vector<Mapped_point> &Qw,
			      int p, const std::vector<double> &u,
			      const std::vector<double> &U,
			      std::vector<Mapped_point> &Pw);
NurbsSurface skinned_surface(const std::vector<Nurbs> &C, int q);
#endif
