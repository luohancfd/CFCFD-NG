## \file geom_test.py
## \ingroup libgeom2
## \brief Exercise the Vector3 and Node3 C++ classes from Python
## \author PJ
## \version 30-Dec-2005

from libgeom2 import *
from math import sqrt, pi

def print_nodes():
    print "There are", len(Node.nodeList), "entries in node_list."
    for node in Node.nodeList:
        print node
    return

def write_vrml_file():
    outfile = open("geom_test_py.wrl", "w");
    outfile.write("#VRML V2.0 utf8\n")
    for n in Node.nodeList:
	outfile.write(n.vrml_str(rgb=vrml_pure_red, radius=0.05, 
                                 offset=0.05, font_size=0.1))
    outfile.close()

print "Begin geom_test..."

print "First, exercise the Vector3 class."
a = Vector(1.1, 2.0, 3.0)
print "initial : a=", a
a += 2.0
print "a+=2 : a=", a
b = a
b += a
print "b=a; b+=a : b=", b
print "vabs(b)=", vabs(b)
print "vabs(b.norm)=", vabs(b.norm()), "; b=", b
print "+a=", +a, "; -b=", -b
a -= 2 # implicit conversion
print "a-=2: a=", a
a -= b; print "a-=b: a=", a
c = a - b;
print "c=a-b: c=", c
c = 2.0 * a + 3.0;
print "c=2.0*a+3.0: c=", c # implicit conversion
d = Vector()
d = 2.0 * a + Vector(3.0);
print "d=Vector(2.0)*a+Vector(3.0): d=", d
print "equal(a,b,1.0e-6)=", equal(a,b,1.0e-6)
print "equal(c,d)=", equal(c,d)
print "equal(c,1.0)=", equal(c,Vector(1.0)) # implicit conversion not available
a = Vector(sqrt(2),sqrt(2),0.0)
b = Vector(-a.y, a.x, 0.0)
print "a=", a, " b=", b, " cross(a,b)=", cross(a,b)
print "unit(cross(a,b))=", unit(cross(a,b))

p = Vector(1.0, 1.0, 0.0)
n = Vector(1.0, 1.0, 0.0)
a.mirror_image(p, n)
print "a.mirror_image(", p, ",", n, ")= ", a

print "Test some of the geometry functions..."
na = Node(a, label="A")
nb = Node(b, label="B")
print "na=", na, " nb=", nb
print_nodes()
write_vrml_file()
print "equal(na,nb)=", equal(na, nb)
na.translate(10.0); nb.translate(d)
print "After translating:"; print_nodes()
na.rotate_about_zaxis(pi);
nb.rotate_about_zaxis(pi/2).rotate_about_zaxis(pi/2)
print "After rotating 180 degrees:"; print_nodes()
nc = Node(na); nc.label = "C"
nd = nb; nd.label = "D" # Note that nd references the original nb object.
ne = na.copy(); ne.label = "E"
print "After copying:"; print_nodes();

print "Projection onto a plane and testing within a triangle."
a = Vector(1.0, 0.0, 0.0); b = Vector(0.0, 1.0, 0.0); c = Vector(0.0, 0.0, 1.0)
q = Vector(0.0, 0.0, 0.0); qr = Vector(0.5, 0.5, 0.5);
result_flag = project_onto_plane(q, qr, a, b, c);
print "q=", q, "result_flag=", result_flag, \
      "inside_triangle=", inside_triangle(q, a, b, c);
q = Vector(0.0, 0.5, 0.5)
print "q=", q, "inside_triangle=", inside_triangle(q, a, b, c);
q = Vector(0.0, 0.51, 0.5)
print "q=", q, "inside_triangle=", inside_triangle(q, a, b, c);

print "Properties of geometric primitives."
p0 = Vector3(0.0, 0.0, 0.0)
p1 = Vector3(1.0, 0.0, 0.0)
p2 = Vector3(1.0, 1.0, 0.0)
p3 = Vector3(0.0, 1.0, 0.0)
p4 = Vector3(0.0, 0.0, 1.0)
p5 = Vector3(1.0, 0.0, 1.0)
p6 = Vector3(1.0, 1.0, 1.0)
p7 = Vector3(0.0, 1.0, 1.0)
volume = hexahedron_volume(p0, p1, p2, p3, p4, p5, p6, p7)
centroid = hexahedron_centroid(p0, p1, p2, p3, p4, p5, p6, p7)
print "hexahedron: centroid=", centroid, ", volume=", volume

area = quad_area(p0, p1, p2, p3)
centroid = quad_centroid(p0, p1, p2, p3)
n = quad_normal(p0, p1, p2, p3)
t1 = quad_tangent1(p0, p1, p2, p3)
t2 = quad_tangent2(p0, p1, p2, p3)
print "quadrilateral: centroid=", centroid, ", area=", area
print "               n=", n, ",t1=", t1, ", t2=", t2

print "Done."

