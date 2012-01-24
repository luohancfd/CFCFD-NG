/** \file surface_test.cxx
 *  \ingroup libgeom2
 *  \brief Exercise the Surface classes.
 *  \author PJ
 *  \version 09-Jan-2006
 *
 */
#include <iostream>
#include <math.h>
#include "geom.hh"
#include "gpath.hh"
#include "surface.hh"
using namespace std;


int main()
{
    cout << "\n\n-----------------------------------------------" << endl;
    cout << "Begin surface_test..." << endl;

    cout << "\nSimple Coons surface:" << endl;
    Line south = Line(Vector3(-1.0, -1.0, 0.0), Vector3(1.0, -0.7, 0.0));
    Line north = Line(Vector3(-1.0, 1.0, 0.0), Vector3(1.0, 1.0, 0.0));
    Line west = Line(Vector3(-1.0, -1.0, 0.0), Vector3(-1.0, 1.0, 0.0));
    Line east = Line(Vector3(1.0, -0.7, 0.0), Vector3(1.0, 1.0, 0.0));
    AOPatch surf1 = AOPatch(south, north, west, east, "AO_SURFACE");
    cout << surf1 << endl;
    cout << "surf1.eval(0.25,0.75)=" << surf1.eval(0.25,0.75) << endl;

    cout << "\nSimple Cylindrical (Coons) surface:" << endl;
    Line c0 = Line(Vector3(1.0, 0.0, 0.0), Vector3(1.0, 0.0, 1.0));
    Line c1 = c0; c1.translate(-1.0, 1.0, 0.0);
    Arc c2 = Arc(Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(0,0,0));
    Arc c3 = c2; c3.translate(0.0, 0.0, 1.0);
    CoonsPatch surf2 = CoonsPatch(c0, c1, c2, c3, "Cylinder_surface");
    cout << surf2 << endl;
    cout << "surf2.eval(0.5,0.5)=" << surf2.eval(0.5,0.5) << endl;

    Vector3 p00 = Vector3(-1.0, -1.0, 2.0);
    Vector3 p10 = Vector3(1.0, -1.0, 2.0);
    Vector3 p11 = Vector3(1.0, 1.0, 2.0);
    Vector3 p01 = Vector3(-1.0, 1.0, 2.0);
    CoonsPatch surf3 = CoonsPatch(p00, p10, p11, p01, "Square_surface");
    cout << surf3 << endl;
    cout << "surf3.eval(0.25,0.75)=" << surf3.eval(0.25,0.75) << endl;

    cout << "\nPath on the surface (reverse diagonal):" << endl;
    PathOnSurface psurf = PathOnSurface(surf3, LinearFunction(), 
					LinearFunction(-1.0,1.0));
    cout << "psurf=" << psurf << endl;
    int n = 5;
    double dt = 1.0/n;
    cout << "    t     position(t)" << endl;
    for ( int i = 0; i <= n; ++i ) {
	double t = dt * i;
	cout << "pos(" << t << ")=" << psurf.eval(t) << endl; 
    }
    cout << "Done." << endl;
    return 0;
}
