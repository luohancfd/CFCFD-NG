## \file gpath_test_bezier.py
## \ingroup libgeom2
## \brief Exercise the Path C++ classes for Bezier curves from Python
## \author PJ
## \version 19-Sep-2012

from libgeom2 import *
from math import sqrt, pi

print "Begin gpath_test_bezier.py -- to test arc_length parameterization ..."

print "Bezier curves (approximating straight lines):"
Bv = [Vector3(1.0, 1.0, 0.0), Vector3(2.0, 2.0, 0.0),
      Vector3(3.0, 3.0, 0.0), Vector3(4.0, 4.0, 0.0)]
bez0 = Bezier(Bv, "BEZIER0")
print "bez0=", bez0
print "bez0.eval(0.5)=", bez0.eval(0.5), " length=", bez0.length()
Bv = [Vector3(1.0, 1.0, 0.0), Vector3(1.5, 1.5, 0.0),
      Vector3(2.0, 2.0, 0.0), Vector3(4.0, 4.0, 0.0)]
bez1 = Bezier(Bv, "BEZIER1", 0.0, 1.0, 1)
print "bez1=", bez1
print "bez1.eval(0.5)=", bez1.eval(0.5), " length=", bez1.length()
bez2 = bez1.copy(direction=-1)
for i in range(5):
    t = i * 0.25
    print "%10f %30s %45s %45s" % (t, bez0.eval(t), bez1.eval(t), bez2.eval(t))

bez0.add_point(Vector3(5.0, 5.0, 0.0)).add_point(Vector3(6.0, 6.0, 0.0))
print "bez0=", bez0
p2 = bez0.eval(0.5)
print "bez0.eval(0.5)=", p2, " length=", bez0.length()
bez1.add_point(Vector3(4.5, 4.5, 0.0)).add_point(Vector3(6.0, 6.0, 0.0))
print "bez1=", bez1
p2 = bez1.eval(0.5)
print "bez1.eval(0.5)=", p2, " length=", bez1.length()


print "Done."

