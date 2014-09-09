/** \file gpath_test.cxx
 *  \ingroup libgeom2
 *  \brief Exercise the geometric-path classes.
 *  \author PJ
 *  \version 27-Dec-2005
 *  \version 17-Jan-2006 remove C++ Node3 class
 *
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include "geom.hh"
#include "gpath.hh"
#include "gpath_utils.hh"
using namespace std;

// We want a vector container (from the standard library) 
// to hold pointers to the Nodes.
static vector<Vector3*> node_list;
static vector<Path*> path_list;

void print_nodes( ostream &os=cout )
{
    int i, n;
    n = node_list.size();
    os << "node_list contains " << n << " entries" << endl;
    for ( i = 0; i < n; ++i ) {
	// print what is pointed to...
	os << *(node_list[i]) << endl;
    }
}

int main()
{
    cout << "\n\n-----------------------------------------------" << endl;
    cout << "Begin gpath_test..." << endl;

    cout << "\nStraight lines:" << endl;
    Vector3 a(1.0,2.0,3.0); node_list.push_back(&a);
    Vector3 b = Vector3(2.0, 3.0, 4.0); node_list.push_back(&b);
    print_nodes();
    cout << "end points : a=" << a << "  b=" << b << endl;
    Line ab = Line(a,b, "AB"); path_list.push_back(&ab);
    cout << "ab=" << ab << endl;
    print_nodes();
    cout << "ab.eval(0.5)=" << ab.eval(0.5) << endl;
    print_nodes();
    ab.reverse();
    cout << "after reversing: ab=" << ab << endl;
    print_nodes();

    cout << "\nArcs" << endl;
    Vector3 a1 = Vector3(0.0, 0.0, 1.0); node_list.push_back(&a1);
    Vector3 b1 = Vector3(1.0, 0.0, 0.0); node_list.push_back(&b1);
    Vector3 c = Vector3(0.0, 0.0, 0.0);
    Arc abc = Arc(a1, b1, c, "ABC");
    path_list.push_back(&abc);
    cout << "abc=" << abc << endl;
    Vector3 p = abc.eval(0.5);
    cout << "abc.eval(0.5)=" << p << " R=" << vabs(p-c) << endl; 
    abc.reverse();
    cout << "after reversing: abc=" << abc << endl;
    print_nodes();

    cout << "\nArc3s" << endl;
    Vector3 start = a1; Vector3 mid = p; Vector3 end = b1;
    Arc3 abc3 = Arc3(start, mid, end, "ABC3");
    cout << "abc3=" << abc3 << endl;
    Vector3 p2 = abc3.eval(0.5);
    cout << "abc3.eval(0.5)=" << p2 
	 << " R=" << vabs(p2-abc3.c)
	 << " length=" << abc3.length() 
	 << endl; 
    abc3.reverse();
    cout << "after reversing: abc3=" << abc3 << endl;
    print_nodes();

    cout << "\nBezier curves (approximating circular arcs):" << endl;
    vector<Vector3> Bv;
    double k = 4.0/3.0*(sqrt(2.0) - 1.0);
    Bv.push_back(Vector3(1.0, 0.0, 0.0));
    Bv.push_back(Vector3(1.0,   k, 0.0));
    Bv.push_back(Vector3(  k, 1.0, 0.0));
    Bv.push_back(Vector3(0.0, 1.0, 0.0));
    Bezier bez = Bezier(Bv, "BEZIER");
    cout << "bez=" << bez << endl;
    p2 = bez.eval(0.5);
    cout << "bez.eval(0.5)=" << p2 << " length=" << bez.length() << endl;
    cout << "\nAdd a few more nodes (as a mirror-image in the y-axis)." << endl;
    cout << "Note that we do not expect this high-order Bezier to be a good approximation."
	 << endl;
    bez.add_point(Vector3( 0.0, 1.0, 0.0));
    bez.add_point(Vector3(  -k, 1.0, 0.0));
    bez.add_point(Vector3(-1.0,   k, 0.0));
    bez.add_point(Vector3(-1.0, 0.0, 0.0));
    cout << "bez=" << bez << endl;
    p2 = bez.eval(0.5);
    cout << "bez.eval(0.5)=" << p2 << " length=" << bez.length() << endl;

    cout << "\nA Nurbs curve (Ex4.1 from Piegl and Tiller, 1997):" << endl;
    vector<Vector3> P;
    P.push_back(Vector3(0.0, 0.0));
    P.push_back(Vector3(1.0, 1.0));
    P.push_back(Vector3(3.0, 2.0));
    P.push_back(Vector3(4.0, 1.0));
    P.push_back(Vector3(5.0, -1.0));
    // All weights 1.0, except w[1] = 4.0
    vector<double> w(5, 1.0);
    w[1] = 4.0;
    int deg = 2;
    vector<double> U(8, 0.0);
    U[0] = 0.0; U[1] = 0.0; U[2] = 0.0; U[3] = 1.0;
    U[4] = 2.0; U[5] = 3.0; U[6] = 3.0; U[7] = 3.0;
    double t = 1./3.; // equivalent to u = 1.0 with this
                      // knot vector
    Nurbs nurbs(P, w, deg, U);
    cout << "nurbs=" << nurbs << endl;
    p2 = nurbs.eval(t);
    cout << "nurbs.eval(1./3.)=" << p2 << endl;
    cout << "Expected value= (7/5, 6/5, 0)\n";
    cout << "Check that the Nurbs curve can correctly represent the Bezier.\n";
    cout << "Based on the Bezier used earlier, approximating a circuar arc.\n";
    // Set all weights to 1.0;
    w.resize(4);
    w[0] = 1.0; w[1] = 1.0; w[2] = 1.0; w[3] = 1.0;
    // Set knot vector for a degree 3 Bezier.
    U[0] = 0.0; U[1] = 0.0; U[2] = 0.0; U[3] = 0.0;
    U[4] = 1.0; U[5] = 1.0; U[6] = 1.0; U[7] = 1.0;
    deg = 3;
    Nurbs nurbs2(Bv, w, deg, U);
    cout << "nurbs2=" << nurbs2 << endl;
    p2 = nurbs2.eval(0.5);
    cout << "nurbs2.eval(0.5)=" << p2 << " length=" << nurbs2.length() << endl;

    cout << "\nPolylines:" << endl;
    vector<Path*> segments;
    Bezier bez1 = Bezier(Bv, "BEZIER"); 
    segments.push_back(&bez1);
    abc = Arc(Vector3(0.0,1.0,0.0), Vector3(-1.0,0.0,0.0), Vector3(0.0,0.0,0.0));
    segments.push_back(&abc); 
    Polyline pline = Polyline(segments, "POLYLINE");
    cout << "pline=" << pline << endl;
    Vector3 p3 = pline.eval(0.5);
    cout << "pline.eval(0.5)=" << p3 << " length=" << pline.length() << endl;

    cout << "\nTranslate the Polyline" << endl;
    pline.translate(-2.0, 0.0, 0.0);
    cout << "pline=" << pline << endl;
    p3 = pline.eval(0.5);
    cout << "pline.eval(0.5)=" << p3 << " length=" << pline.length() << endl;

    cout << "\nPolyline within a Polyline" << endl;
    segments.push_back(&pline);
    Polyline pline2 = Polyline(segments, "POLYLINE2");
    cout << "pline2=" << pline2 << endl;
    cout << "pline2.length()=" << pline2.length() << endl;
    int n = 8;
    double dt = 1.0 / n;
    for ( int i = 0; i <= n; ++i ) {
	double t = i * dt;
	cout << "pline2.eval(" << t << ")=" << pline2.eval(t) << endl;
    }

    cout << "\nPolyline within a Polyline -- try subrange." << endl;
    pline2.t0 = 0.25; pline2.t1 = 0.75;
    cout << "pline2=" << pline2 << endl;
    cout << "pline2.length()=" << pline2.length() << endl;
    for ( int i = 0; i <= n; ++i ) {
	double t = i * dt;
	cout << "pline2.eval(" << t << ")=" << pline2.eval(t) << endl;
    }

    cout << "\nSpline" << endl;
    vector<Vector3> ip(5);
    k = 1.0/sqrt(2.0);
    ip[0] = Vector3(1.0, 0.0, 0.0);
    ip[1] = Vector3(k, k, 0.0);
    ip[2] = Vector3(0.0, 1.0, 0.0);
    ip[3] = Vector3(-k, k, 0.0);
    ip[4] = Vector3(-1.0, 0.0, 0.0);
    Spline spl = Spline(ip, "Spline of circle.");
    cout << "spl=" << spl << " spl.length()=" << spl.length() << endl;
    for ( int i = 0; i <= n; ++i ) {
	double t = i * dt;
	p3 = spl.eval(t);
	cout << "spl.eval(" << t << ")=" << p3
	     << " R=" << vabs(p3)
	     << endl;
    }

    cout << "\nMirror image of Spline in x-axis" << endl;
    spl.mirror_image(Vector3(0.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0));
    for ( int i = 0; i <= n; ++i ) {
	double t = i * dt;
	p3 = spl.eval(t);
	cout << "spl.eval(" << t << ")=" << p3
	     << " R=" << vabs(p3)
	     << endl;
    }

    cout << "Done." << endl;

    cout << "\nTesting dist_point_projection routine." << endl;
    Line l = Line(Vector3(0.0, 0.0), Vector3(1.0, 1.0));
    p = Vector3(1.0, 0.0);
    double t_found;
    Vector3 C_found;
    double dist = dist_point_projection(p, l, t_found, C_found);
    cout << "Test point is: " << p << endl;
    cout << "Closest distance to path is: " << dist << endl;
    cout << "Parameter and point on path is: " << t_found << " " << C_found << endl;

    return 0;
}
