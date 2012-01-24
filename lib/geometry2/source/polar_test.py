## \file polar_test.py
## \ingroup libgeom2
## \brief Exercise the PolarPath and PolarSurface for Hannes and Paul.
## \author PJ
## \version 05-Feb-2008

from libgeom2 import *
from math import sqrt

print "---------------------------------------------------"
print "Begin polar_test.py ..."

# Corners for a square in a vertical plane.
a = Vector(0.5, 0.0, 1.0)
dz = Vector(0.0, 0.0, 1.0)
b = a + dz
dy = Vector(0.0, 1.0, 0.0)
c = b + dy
d = c - dz

H = 1.5
ab = PolarPath(Line(a,b), H)
bc = PolarPath(Line(b,c), H)
cd = PolarPath(Line(c,d), H)
da = PolarPath(Line(d,a), H)

p = ab.eval(0.5);
print "ab=", ab, "ab.eval(0.5)=", p

# Now test surface construction.
a = Vector(0.5, -1.0, H)
dy = Vector(0.0, 2.0, 0.0)
b = a + dy
dx = Vector(1.0, 0.0, 0.0)
c = b + dx
d = c - dy
abcd = PolarSurface(CoonsPatch(a, b, c, d), H)

p = abcd.eval(0.5, 0.5);
print "abcd=", abcd, "abcd.eval(0.5,0.5)=", p

print "Render to VRML"
outfile = open("polar_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
for pth in [ab, bc, cd, da]:
    outfile.write(pth.vrml_str() + "\n")
outfile.write(abcd.vrml_str() + "\n")
outfile.close()

print "Done."

