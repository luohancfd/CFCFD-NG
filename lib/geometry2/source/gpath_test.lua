-- This test file mimics the tests in:
--   gpath_test.cxx
--   gpath_test.py
--
-- This source is based on the code by
-- P. A. Jacobs in gpath_test.cxx and gpath_test.py
--
-- Author: Rowan J. Gollan
-- Date: 29-Mar-2008
--

print("Begin gpath_test...\n")
print("\nStraight lines:")
a = Vector3(1.0, 2.0, 3.0)
b = Vector3(2.0, 3.0, 4.0)
print("end points : a=", a, " b=", b)
ab = Line(a, b)
print("ab=", ab)
print("ab:eval(0.5)=", ab:eval(0.5))
ab:reverse()
print("After reversing: ab=", ab)

print("\nArcs:")
a1 = Vector3(0.0, 0.0, 1.0)
b1 = Vector3(1.0, 0.0, 0.0)
c = Vector3(0.0, 0.0, 0.0)
abc = Arc(a1, b1, c, {label="ABC"})
print("abc=", abc)
p = abc:eval(0.5)
print("abc:eval(0.5)=", p, " R=", vabs(p-c))
abc:reverse()
print("after reversing: abc=", abc)

print("\nArc3s:")
start = a1; mid = p; finish = b1
abc3 = Arc3(start, mid, finish, {label="ABC3"})
print("abc3=", abc3)
p2 = abc3:eval(0.5)
print("abc3:eval(0.5)", p2,
      " R=", vabs(p2-abc3.c),
      " length=", abc3:length())
abc3:reverse()
print("after reversing: abc3=", abc3)

print("\nBezier curves (approximating circular arcs):")
Bv = {}
k = 4.0/3.0*(math.sqrt(2.0) - 1.0)
Bv[#Bv + 1] = Vector3(1.0, 0.0, 0.0)
Bv[#Bv + 1] = Vector3(1.0,   k, 0.0)
Bv[#Bv + 1] = Vector3(  k, 1.0, 0.0)
Bv[#Bv + 1] = Vector3(0.0, 1.0, 0.0)
bez = Bezier(Bv, {label="BEZIER"})
print("bez=", bez)
p2 = bez:eval(0.5)
print("bez:eval(0.5)=", p2, " length=", bez:length())
print("\nAdd a few more nodes (as a mirror image in the y-axis).")
print("Note that we do not expect this high-order Bezier to be a good approximation.")
bez:add_point(Vector3( 0.0, 1.0, 0.0))
bez:add_point(Vector3(  -k, 1.0, 0.0))
bez:add_point(Vector3(-1.0,   k, 0.0))
bez:add_point(Vector3(-1.0, 0.0, 0.0))
print("bez=", bez)
p2 = bez:eval(0.5)
print("bez:eval(0.5)=", p2, " length=", bez:length())
print("bez:B(1)=", bez:B(1))

print("\nPolylines:")
segments = {}
bez1 = Bezier(Bv, {label="BEZIER"})
segments[#segments + 1] = bez1
abc = Arc(Vector3(0.0, 1.0, 0.0), Vector3(-1.0, 0.0, 0.0), Vector3(0.0, 0.0, 0.0))
segments[#segments + 1] = abc
pline = Polyline(segments, {label="POLYLINE"})
print("pline=", pline)
p3 = pline:eval(0.5)
print("pline:eval(0.5)=", p3, " length=", pline:length())

print("\nTranslate to Polyline")
pline:translate(-2.0, 0.0, 0.0)
print("pline=", pline)
p3 = pline:eval(0.5)
print("pline:eval(0.5)=", p3, " length=", pline:length())

print("\nPolyline within a Polyline")
segments[#segments + 1] = pline
pline2 = Polyline(segments, "POLYLINE2")
print("pline2=", pline2)
print("pline2:length()=", pline2:length())
n = 8
dt = 1/n
for t=0,1,dt do
   print("pline2:eval(", t, ")=", pline2:eval(t))
end

print("\nPolyline within a Polyline --- try a subrange.")
pline2.t0 = 0.25; pline2.t1 = 0.75
print("pline2=", pline2)
print("pline2:length()=", pline2:length())
for t=0,1,dt do
   print("pline2:eval(", t, ")=", pline2:eval(t))
end

print("\nSpline:")
ip = {}
k = 1.0/math.sqrt(2.0)
ip[1] = Vector3(1.0, 0.0, 0.0)
ip[2] = Vector3(k, k, 0.0)
ip[3] = Vector3(0.0, 1.0, 0.0)
ip[4] = Vector3(-k, k, 0.0)
ip[5] = Vector3(-1.0, 0.0, 0.0)
spl = Spline(ip, {label="Spline of circle."})
print("spl=", spl, " spl:length()=", spl:length())
for t=0,1,dt do
   p3 = spl:eval(t)
   print("spl:eval(", t, ")=", p3, " R=", vabs(p3))
end

print("\nMirror image of Spline in x-axis")
spl:mirror_image(Vector3(0.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0))
for t=0,1,dt do
   p3 = spl:eval(t)
   print("spl:eval(", t, ")=", p3, " R=", vabs(p3))
end

print("\nXPoly:")
print("Simple line: y = 0.0 + 1.0*x")
print("x0 = 0.0, x1 = 1.0, B = {0.0, 1.0}")
a = XPoly(0.0, 1.0, {0.0, 1.0})
print("xpoly=", a)
print("XPoly:xeval(0.5): expect Vector3(0.5, 0.5, 0.0)")
print(a:xeval(0.5))
print("For this simple line, eval should be equivalent to xeval.")
print("XPoly:eval(0.5): expect Vector3(0.5, 0.5, 0.0)")
print(a:eval(0.5))
print("XPoly:dydx(0.2): expect 1.0 --- the gradient")
print(a:dydx(0.2))
print("XPoly:dpdt(0.4): expect Vector3(1.0, 1.0, 0.0)")
print(a:dpdt(0.4))

print("\nNow try a parabola: y = 2.0 + 3.0*x + 4.0*x^2")
print("x0 = -5.0, x1 = 10.0, B = {2.0, 3.0, 4.0}")
b = XPoly(-5.0, 10.0, {2.0, 3.0, 4.0})
print("XPoly:xeval(4.0): expect Vector3(4, 78, 0)")
print(b:xeval(4.0))
print("XPoly:eval(0.6): expect Vector3(4, 78, 0)")
print(b:eval(0.6))
print("XPoly:dydx(-5.0): expect -37")
print(b:dydx(-5.0))
print("XPoly:dpdt(1.0): expect Vector3(15, 1245, 0)")
print(b:dpdt(1.0))

print("\nXBezier:")
print("Test printing of XBezier.")
xbez = XBezier({Vector3(0.0, 0.0, 0.0),
		Vector3(1.0, 2.0, 0.0),
		Vector3(2.0, 4.0, 0.0),
		Vector3(3.0, 6.0, 0.0)})
print("xbez=", xbez)
print("xbez:xeval(3.0)= ", xbez:xeval(3.0), " xbez:length()=", xbez:length())

print("Done.") 
