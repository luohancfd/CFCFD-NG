/** \file revolved_surface_test.cxx
 *  \ingroup libgeom2
 *  \brief Exercise the RevolvedSurface and MappedSurface classes.
 *  \author PJ
 *  \version 28-Jan-2006
 *
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "geom.hh"
#include "gpath.hh"
#include "surface.hh"
using namespace std;


int main()
{
    cout << "\n\n-----------------------------------------------" << endl;
    cout << "Begin revolved_surface_test..." << endl;

    cout << "\nSimple RevolvedSurface:" << endl;
    Line path = Line(Vector3(0.0,0.0,0.0), Vector3(1.0,1.0, 0.0));
    RevolvedSurface surf1 = RevolvedSurface(&path, "REVOLVED_SURFACE");
    cout << "surf1=" << surf1 << endl;
    cout << "surf1.eval(0.25,0.75)=" << surf1.eval(0.25,0.75) << endl;

    cout << "\nMapped surface:" << endl;
    Vector3 p0(-0.5, -0.5, -0.5);
    Vector3 p1(-0.5, -0.5, 0.5);
    Vector3 p2(-0.5, 0.5, 0.5);
    Vector3 p3(-0.5, 0.5, -0.5);
    CoonsPatch qsurf = CoonsPatch(p0, p1, p2, p3, "QUERY_SURFACE");
    MappedSurface msurf = MappedSurface(&qsurf, &surf1);
    cout << "msurf=" << msurf << endl;
    cout << "msurf.eval(0.25,0.75)=" << msurf.eval(0.25,0.75) << endl;

    return 0;
}
