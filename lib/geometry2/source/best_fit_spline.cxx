/** \file spline_fitter.cxx
 *
 *  \author Rowan J. Gollan
 *  \version 15-Jan-2007
 *
 **/
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

#include "../../nm/source/fobject.hh"
#include "geom.hh"
#include "gpath.hh"
#include "best_fit_spline.hh"

using namespace std;

Best_fit_spline::
Best_fit_spline(string input_file,
		double errwt1, double errwt2,
		double y0, double x1)
    : MultivariateFunction(), errwt1_(errwt1), errwt2_(errwt2),
      y0_(y0), x1_(x1)
{
    double x, y;
    ifstream infile;
    infile.open(input_file.c_str());
    if( infile.fail() ) {
	cout << "Could not open " << input_file << "; Bailing out!\n";
	exit(1);
    }

    while( ! infile.eof() ) {
	infile >> x;
	infile >> y;
	if( infile.eof() )
	    break;
	if( x > x1_ ) continue;
	points_.push_back(Vector3(x, y));
    }

    infile.close();
}

Best_fit_spline::
Best_fit_spline(const Best_fit_spline &bs)
    : errwt1_(bs.errwt1_), errwt2_(bs.errwt2_),
      y0_(bs.y0_), x1_(bs.x1_), points_(bs.points_) {}
 

Best_fit_spline::
~Best_fit_spline() {}

Best_fit_spline*
Best_fit_spline::
clone() const
{
    return new Best_fit_spline(*this);
}

double 
Best_fit_spline::
eval( vector<double> &x )
{
    vector<Vector3> p;
    p.push_back(Vector3(x[0], y0_));
    for(size_t i = 1; i < x.size()-1; i += 2) {
	p.push_back(Vector3(x[i],x[i+1]));
    }
    p.push_back(Vector3(x1_, x[x.size()-1]));

    Spline mySpline(p);

    double err1 = 0.0;
    double y_val, x_val;
    Vector3 test;
    for(size_t i = 0; i < points_.size()-1; ++i ) {
	x_val = points_[i].x;
	y_val = points_[i].y;
	if( y_val < y0_ ) continue;
	test = mySpline.eval_from_y(y_val);
	err1 += (x_val - test.x)*(x_val - test.x);
    }
    // Additional error if end point is straying.
    double end_s_y = mySpline.eval(1.0).y;
    double end_y = points_[points_.size()-1].y;
    err1 += (end_y - end_s_y)*(end_y - end_s_y);


//     int nsp = p.size();
//     double err2 = 0.0;
//     double ideal_length = mySpline.length() / (nsp-1);
//     vector<Vector3> ideal_points;
//     Path *myPath = &mySpline;
//     for(int i = 0; i < (nsp-2); ++i ) {
// 	ideal_points.push_back(myPath->point_from_length((i+1)*ideal_length));
//     }

//     Vector3 actual_point;
//     for(int i = 0; i < (nsp-2); ++i ) {
// 	actual_point = mySpline.eval(mySpline.t_seg[i]);
// 	cout << "i= " << i << " t= " << mySpline.t_seg[i] << " ideal= " << ideal_points[i] << " actual= " << actual_point << endl;
// 	err2 += vabs(ideal_points[i] - actual_point);
// 	cout << "error= " << vabs(ideal_points[i] - actual_point) << endl;
//     }

    double error = err1;
    //errwt1_*err1 + errwt2_*err2;

    return error;
}

string
Best_fit_spline::
str() const
{
    return string("Best_fit_spline() object.");
}

