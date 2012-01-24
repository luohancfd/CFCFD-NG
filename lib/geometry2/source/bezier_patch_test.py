## \file bezier_patch_test.py
## \ingroup libgeom2
## \brief Exercise the BezierPatch surface class.
## \author PJ
## \version 07-Mar-2006

from math import cos, sqrt
from libgeom2 import *

print "\n\n-----------------------------------------------"
print "Begin bezier_patch_test..."

print "\nSimple BezierPatch surface:"
origin = Vector(0.0, 0.0, 0.0)
dx = 0.1; dy = 0.2
n = 6; m = 5
Q = []
for i in range(n+1):
    for j in range(m+1):
        x = i * dx
        y = j * dy
        z = cos(sqrt(x*x + y*y))
        Q.append(origin + Vector(i*dx, j*dy, z))
surf1 = BezierPatch(Q, n, m, "BEZIER_PATCH_SURFACE")
print surf1
p = surf1.eval(0.25,0.75)
print "surf1.eval(0.25,0.75)=", p
print "expect z=", cos(sqrt(p.x*p.x  +p.y*p.y)), "very roughly"

print "Render to VRML"
outfile = open("bezier_patch_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(surf1.vrml_str() + "\n")
outfile.close()

print "Done."
