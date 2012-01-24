/** \file triangle_patch_test.cxx
 *  \ingroup libgeom2
 *  \brief Exercise the TrianglePatch class.
 *  \author PJ
 *  \version 26-Jan-2006
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
    cout << "Begin triangle_patch_test..." << endl;

    cout << "\nSimple TrianglePatch surface:" << endl;
    vector<Vector3*> p = vector<Vector3*>();
    p.push_back(new Vector3(0.0,0.0,0.0));
    p.push_back(new Vector3(1.0,0.0,0.0));
    p.push_back(new Vector3(1.0,1.0,0.0));
    p.push_back(new Vector3(0.0,1.0,0.0));
    p.push_back(new Vector3(0.5,0.5,1.0));
    vector<int>itri = vector<int>();
    itri.push_back(0); itri.push_back(1); itri.push_back(4);
    itri.push_back(1); itri.push_back(2); itri.push_back(4);
    itri.push_back(2); itri.push_back(3); itri.push_back(4);
    itri.push_back(3); itri.push_back(0); itri.push_back(4);
    vector<int>iA = vector<int>();
    iA.push_back(0); iA.push_back(1);
    vector<int>iB = vector<int>();
    iB.push_back(3); iB.push_back(2);
    vector<int>iC = vector<int>();
    iC.push_back(0); iC.push_back(3);
    vector<int>iD = vector<int>();
    iD.push_back(1); iD.push_back(2);
    TrianglePatch surf1 = TrianglePatch(p, itri, iA, iB, iC, iD, "TRI_SURFACE");
    cout << "surf1=" << surf1 << endl;
    cout << "surf1.eval(0.25,0.75)=" << surf1.eval(0.25,0.75) << endl;

    cout << "\nSimple Cylindrical (Coons) surface:" << endl;
    Line c0 = Line(Vector3(1.0, 0.0, 0.0), Vector3(1.0, 0.0, 1.0));
    Line c1 = c0; c1.translate(-1.0, 1.0, 0.0);
    Arc c2 = Arc(Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(0,0,0));
    Arc c3 = c2; c3.translate(0.0, 0.0, 1.0);
    CoonsPatch surf2 = CoonsPatch(c0, c1, c2, c3, "Cylinder_surface");
    cout << surf2 << endl;
    cout << "surf2.eval(0.5,0.5)=" << surf2.eval(0.5,0.5) << endl;

    cout << "\nTriangulated Cylindrical surface:" << endl;
    TrianglePatch surf3 = TrianglePatch( &surf2, 1, 5, "triangulated-patch");
    cout << surf3 << endl;
    cout << "surf3.eval(0.5,0.5)=" << surf3.eval(0.5,0.5) << endl;

    // Clean up allocated memory so that valgrind doesn't complain.
    for ( size_t i = 0; i < p.size(); ++i ) delete p[i];
    cout << "Done." << endl;
    return 0;
}
