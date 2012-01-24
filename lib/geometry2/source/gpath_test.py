## \file gpath_test.py
## \ingroup libgeom2
## \brief Exercise the Path C++ classes from Python
## \author PJ
## \version 01-Jan-2006

from libgeom2 import *
from math import sqrt, pi

def print_nodes():
    print "There are", len(Node.nodeList), "entries in node_list."
    for node in Node.nodeList:
        print node
    return

print "---------------------------------------------------"
print "Begin gpath_test.py ..."

print "\nStraight lines:"
a = Node(1.0,2.0,3.0, "A")
b = Node(2.0, 3.0, 4.0, "B")
print_nodes()
print "end points : a=", a, "  b=", b
ab = Line(a,b, "AB");
print "ab=", ab
print_nodes()
print "ab.eval(0.5)=", ab.eval(0.5)
print_nodes()
ab.reverse()
print "after reversing: ab=", ab
print_nodes()
ab.rotate_about_zaxis(pi/2).rotate_about_zaxis(pi/2)
print "after rotating about the zaxis 180degrees: ab=", ab
print_nodes()

print "Render to VRML"
outfile = open("gpath_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(a.vrml_str(radius=0.05) + "\n")
outfile.write(b.vrml_str(radius=0.05) + "\n")
outfile.write(ab.vrml_str() + "\n")

print "\nArcs:"
a = Node(0.0, 0.0, 1.0, "a") 
b = Node(1.0, 0.0, 0.0, "b")
c = Vector3(0.0, 0.0, 0.0)
abc = Arc(a, b, c, "ABC")
print "abc=", abc
p = abc.eval(0.5);
print "abc.eval(0.5)=", p, \
      " R=", vabs(p-c), \
      " length=", abc.length()
print "abc.dpdt(0.5)=", abc.dpdt(0.5), \
      "abc.dpdt(0.0)=", abc.dpdt(0.0), \
      "abc.dpdt(1.0)=", abc.dpdt(1.0)
abc.reverse()
print "after reversing: abc=", abc
print_nodes()

print "\nArc3s:"
start = a; mid = p; end = b
abc3 = Arc3(start, mid, end, "ABC3")
print "abc3=", abc3
p2 = abc3.eval(0.5)
print "abc3.eval(0.5)=", p2, \
      " R=", vabs(p2-abc3.c), \
      " length=", abc.length()
abc3.reverse()
print "after reversing: abc3=", abc3
print_nodes()

# Continue the VRML rendering of some elements.
outfile.write(a.vrml_str(radius=0.05) + "\n")
outfile.write(b.vrml_str(radius=0.05) + "\n")
outfile.write(abc.vrml_str() + "\n")
outfile.close()

print "\nBezier curves (approximating circular arcs):"
k = 4.0/3.0*(sqrt(2.0) - 1.0)
Bv = [Vector3(1.0, 0.0, 0.0), Vector3(1.0,   k, 0.0),
      Vector3(  k, 1.0, 0.0), Vector3(0.0, 1.0, 0.0)]
bez = Bezier(Bv, "BEZIER");
print "bez=", bez
p2 = bez.eval(0.5);
print "bez.eval(0.5)=", p2, " length=", bez.length()
print "\nAdd a few more nodes (as a mirror-image in the y-axis)."
print "Note that we do not expect this high-order Bezier to be a good approximation."
bez.add_point(Vector3( 0.0, 1.0, 0.0)).add_point(Vector3(  -k, 1.0, 0.0))
bez.add_point(Vector3(-1.0,   k, 0.0)).add_point(Vector3(-1.0, 0.0, 0.0))
print "bez=", bez
p2 = bez.eval(0.5);
print "bez.eval(0.5)=", p2, " length=", bez.length()

print "\nPolylines:"
bez1 = Bezier(Bv, "BEZIER")
abc = Arc(Vector3(0.0,1.0,0.0), Vector3(-1.0,0.0,0.0), Vector3(0.0,0.0,0.0))
pline = Polyline([bez1, abc], "POLYLINE")
print "pline=", pline
p3 = pline.eval(0.5)
print "pline.eval(0.5)=", p3, " length=", pline.length()

print "\nTranslate the Polyline"
pline.translate(-2.0, 0.0, 0.0)
print "pline=", pline
p3 = pline.eval(0.5)
print "pline.eval(0.5)=", p3, " length=", pline.length()

print "\nPolyline within a Polyline"
pline2 = Polyline([bez1, abc, pline], "POLYLINE2")
print "pline2=", pline2
print "pline2.length=", pline2.length()
n = 8
dt = 1.0 / n
for i in range(n+1):
    t = i * dt
    print "pline2.eval(", t, ")=", pline2.eval(t)

print "\nPolyline within a Polyline -- try subrange"
pline2.t0 = 0.25; pline2.t1 = 0.75
print "pline2=", pline2
print "pline2.length=", pline2.length()
for i in range(n+1):
    t = i * dt
    print "pline2.eval(", t, ")=", pline2.eval(t)

print "\nSpline"
k = 1.0/sqrt(2.0)
ip = [Vector3(1.0, 0.0, 0.0), Vector3(k, k, 0.0), Vector3(0.0, 1.0, 0.0),
      Vector3(-k, k, 0.0), Vector3(-1.0, 0.0, 0.0)]
spl = Spline(ip, "Spline of circle.")
print "spl=", spl, " spl.length()=", spl.length()
for i in range(n+1):
    t = i * dt
    p3 = spl.eval(t)
    print "spl.eval(", t, ")=", p3, " R=", vabs(p3)

print "\nMirror image of Spline in x-axis"
spl.mirror_image(Vector(0.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0))
for i in range(n+1):
    t = i * dt
    p3 = spl.eval(t)
    print "spl.eval(", t, ")=", p3, " R=", vabs(p3)


print "\nXPoly:"
print "Simple line: y = 0.0 + 1.0*x"
print "x0 = 0.0, x1 = 1.0, B = vectord([0.0, 1.0])"
a = XPoly(0.0, 1.0, vectord([0.0, 1.0]))
print "XPoly::xeval(0.5): expect Vector3(0.5, 0.5, 0.0)"
print a.xeval(0.5)
print "For this simple line, eval should be equivalent to xeval."
print "XPoly::eval(0.5): expect Vector3(0.5, 0.5, 0.0)"
print a.eval(0.5)
print "XPoly::dydx(0.2): expect 1.0 --- the gradient"
print a.dydx(0.2)
print "XPoly::dpdt(0.4): expect Vector3(1.0, 1.0, 0.0)"
print a.dpdt(0.4)

print ""
print "Now try a parabola: y = 2.0 + 3.0*x + 4.0*x**2"
print "x0 = -5.0, x1 = 10.0, B = vectord([2.0, 3.0, 4.0])"
b = XPoly(-5.0, 10.0, vectord([2.0, 3.0, 4.0]))
print "XPoly::xeval(4.0): expect Vector3(4.0, 78.0, 0.0)"
print b.xeval(4.0)
print "Find t for given x, x=4.0: XPoly::map_x2t(4.0)"
print b.map_x2t(4.0)
print "XPoly::eval(0.6): expect Vector3(4.0, 78.0, 0.0)"
print b.eval(0.6)
print "XPoly::dydx(-5.0): expect -37.0"
print b.dydx(-5.0)
print "XPoly::dpdt(1.0): expect Vector3(15.0, 1245.0, 0.0)"
print b.dpdt(1.0)

print "\nXBezier"
xbez = XBezier(0.0, 3.0, [Vector3(0.0, 0.0, 0.0),
                          Vector3(1.0, 2.0, 0.0),
                          Vector3(2.0, 4.0, 0.0),
                          Vector3(3.0, 6.0, 0.0)])
print "xbez=", xbez
print "XBezier::xeval(3.0): expect Vector3(3.0, 6.0, 0.0)"
print xbez.xeval(3.0)
print "XBezier::length(): expect 6.7082"
print xbez.length()


print "\nXPolyline:"
print "Try connecting two lines."
print "x0 = 0.0, x1 = 1.0, B = vectord([0.0, 1.0])"
a = XPoly(0.0, 1.0, vectord([0.0, 1.0]))
print "x0 = 1.0, x1 = 4.0, B = vectord([1.0, 2.0])"
b = XPoly(1.0, 4.0, vectord([1.0, 2.0]))
segments = vectorXPathPtr(2)
segments[0] = a
segments[1] = b
c = XPolyline(segments)
print "XPolyline::xeval(0.5): expect Vector3(0.5, 0.5, 0.0)"
print c.xeval(0.5)
print "Xpolyline::dydx(0.5): expect 1.0"
print c.dydx(0.5)
print "XPolyline::xeval(1.0): expect Vector3(1.0, 1.0, 0.0)"
print c.xeval(1.0)
print "Xpolyline::dydx(1.0): expect 1.0"
print c.dydx(1.0)
print "XPolyline::xeval(2.0): expect Vector3(2.0, 5.0, 0.0)"
print c.xeval(2.0) 
print "Xpolyline::dydx(2.0): expect 2.0"
print c.dydx(2.0)


print "\nPyFunctionPath:"
def myfun(t): return (1.0*t,2.0*t,1.0)
pp = PyFunctionPath(myfun, "Straight-line")
print pp
print "PyFunctionPath::eval(0.5): expect Vector3(0.5,1.0,1.0)"
print pp.eval(0.5)

print "Done."

