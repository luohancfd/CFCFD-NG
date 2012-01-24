## \file revolved_surface_test.py
## \ingroup libgeom2
## \brief Exercise the RevolvedSurface class.
## \author PJ
## \version 28-Jan-2006

from math import *
from libgeom2 import *

print "\n\n-----------------------------------------------"
print "Begin revolved_surface_test..."

print "\nSimple RevolvedSurface describing a sphere-cone:"
Rnose = 1.0
Angle = 45.0 * pi / 180.0
Length = 1.0
c = Vector(0.0, 0.0, 0.0)    # centre of radius
a = Vector(-Rnose, 0.0, 0.0)  # tip of nose
b = Vector(-Rnose*cos(Angle), Rnose*sin(Angle), 0.0) # join between sphere and cone
d = b + Length * Vector(cos(Angle), sin(Angle), 0.0) # skirt of cone
path = Polyline([Arc(a,b,c), Line(b,d)])
surf1 = RevolvedSurface(path, "REVOLVED_SURFACE")
print surf1
print "surf1.eval(0.25,0.75)=", surf1.eval(0.25,0.75)

print "\nMapped surface:"
p0 = Vector3(-0.5, -0.5, -0.5)
p1 = Vector3(-0.5, -0.5, 0.5)
p2 = Vector3(-0.5, 0.5, 0.5)
p3 = Vector3(-0.5, 0.5, -0.5)
qsurf = CoonsPatch(p0, p1, p2, p3, "QUERY_SURFACE")
msurf = MappedSurface(qsurf, surf1)
print "msurf=", msurf
print "msurf.eval(0.25,0.75)=", msurf.eval(0.25,0.75)

print "Render to VRML"
outfile = open("revolved_surface_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(surf1.vrml_str() + "\n")
outfile.write(msurf.vrml_str() + "\n")
outfile.close()

print "Done."
