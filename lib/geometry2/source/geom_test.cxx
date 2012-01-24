/** \file geom_test.cxx
 *  \ingroup libgeom2
 *  \brief Exercise the Vector3 class.
 *  \author PJ
 *  \version 27-Dec-2005
 *
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include "geom.hh"
using namespace std;

// We want a vector container (from the standard library) 
// to hold pointers to the Nodes.
static vector<Vector3*> node_list;

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
    cout << "Begin geom_test..." << endl;

    cout << "First, exercise the Vector3 class." << endl;
    Vector3 a(1.1,2.0,3.0);
    cout << "initial : a=" << a << endl;
    a += 2;
    cout << "a+=2 : a=" << a << endl;
    Vector3 b = a;
    b += a;
    cout << "b=a; b+=a : b=" << b << endl;
    cout << "vabs(b)=" << vabs(b) << endl;
    cout << "vabs(b.norm)=" << vabs(b.norm()) << "; b=" << b << endl;
    cout << "+a=" << +a << "; -b=" << -b << endl;
    cout << "a-=2: a=" << (a -= 2) << endl; // implicit conversion
    cout << "a-=b: a=" << (a -= b) << endl;
    Vector3 c = a - b;
    cout << "c=a-b: c=" << c << endl;
    c = 2.0 * a + 3.0;
    cout << "c=2.0*a+3.0: c=" << c << endl; // implicit conversion
    Vector3 d;
    d = 2.0 * a + Vector3(3.0);
    cout << "d=Vector3(2.0)*a+Vector3(3.0): d=" << d << endl;
    cout << "equal(a,b,1.0e-6)=" << equal(a,b,1.0e-6) << endl;
    cout << "equal(c,d)=" << equal(c,d) << endl;
    cout << "equal(c,1.0)=" << equal(c,1.0) << endl; // implicit conversion
    a = Vector3(sqrt(2.0),sqrt(2.0),0.0);
    b = Vector3(-a.y, a.x, 0.0);
    cout << "a=" << a << " b=" << b << " cross(a,b)=" << cross(a,b) << endl;
    cout << "unit(cross(a,b))=" << unit(cross(a,b)) << endl;

    Vector3 p = Vector3(1.0, 1.0, 0.0);
    Vector3 n = Vector3(1.0, 1.0, 0.0);
    a.mirror_image(p, n);
    cout << "a.mirror_image(" << p << ", " << n << ")= " << a << endl;

    cout << "Projection onto a plane and testing within a triangle." << endl;
    a = Vector3(1.0, 0.0, 0.0); b = Vector3(0.0, 1.0, 0.0); c = Vector3(0.0, 0.0, 1.0);
    Vector3 q = Vector3(0.0, 0.0, 0.0); Vector3 qr = Vector3(0.5, 0.5, 0.5);
    int result_flag = project_onto_plane(q, qr, a, b, c);
    cout << "q= " << q << "  result_flag= " << result_flag << endl;
    cout << "inside_tringle= " << inside_triangle(q, a, b, c) << endl;

    cout << "Properties of geometric primitives." << endl;
    Vector3 p0 = Vector3(0.0, 0.0, 0.0);
    Vector3 p1 = Vector3(1.0, 0.0, 0.0);
    Vector3 p2 = Vector3(1.0, 1.0, 0.0);
    Vector3 p3 = Vector3(0.0, 1.0, 0.0);
    Vector3 p4 = Vector3(0.0, 0.0, 1.0);
    Vector3 p5 = Vector3(1.0, 0.0, 1.0);
    Vector3 p6 = Vector3(1.0, 1.0, 1.0);
    Vector3 p7 = Vector3(0.0, 1.0, 1.0);
    Vector3 centroid = Vector3();
    double volume = 0.0;
    result_flag = hexahedron_properties( p0, p1, p2, p3, p4, p5, p6, p7,
					 centroid, volume );
    cout << "hexahedron: result_flag= " << result_flag 
	 << ", centroid= " << centroid 
	 << ", volume= " << volume << endl;

    Vector3 t1, t2;
    double area;
    result_flag = quad_properties( p0, p1, p2, p3, 
				   centroid, n, t1, t2, area );
    cout << "quadrilateral: result_flag= " << result_flag
	 << ", centroid= " << centroid << ", area= " << area << endl;
    cout << "               n= " << n << ",t1= " << t1 << ", t2= " << t2 << endl;

    cout << "Done." << endl;
    return 0;
}
