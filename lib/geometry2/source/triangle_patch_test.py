## \file triangle_patch_test.py
## \ingroup libgeom2
## \brief Exercise the TrianglePatch Surface class.
## \author PJ
## \version 26-Jan-2006

from libgeom2 import *

print "\n\n-----------------------------------------------"
print "Begin triangle_patch_test..."

print "\nSimple TrianglePatch surface:"
p = [Vector(0.0,0.0,0.0), Vector(1.0,0.0,0.0), Vector(1.0,1.0,0.0),
     Vector(0.0,1.0,0.0), Vector(0.5,0.5,1.0)]
surf1 = TrianglePatch(p, [0,1,4,1,2,4,2,3,4,3,0,4],
                      [0,1], [3,2], [0,3], [1,2],
                      "TRI_SURFACE")
print surf1
print "surf1.eval(0.25,0.75)=", surf1.eval(0.25,0.75)

print "\nSimple cylindrical (Coons) surface:"
c0 = Line(Vector3(1.0, 0.0, 0.0), Vector3(1.0, 0.0, 1.0))
c1 = c0.copy(); c1.translate(-1.0, 1.0, 0.0)
c2 = Arc(Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(0,0,0))
c3 = c2.copy(); c3.translate(0.0, 0.0, 1.0)
surf2 = CoonsPatch(c0, c1, c2, c3, "Cylinder_surface")
print surf2
print "surf2.eval(0.5,0.5)=", surf2.eval(0.5,0.5)

print  "\nTriangulated cylindrical surface:"
surf3 = TrianglePatch( surf2, 1, 5, "triangulated-patch")
print surf3
print "surf3.eval(0.5,0.5)=", surf3.eval(0.5,0.5)

print "\nJoined triangulated surface."
surf2.translate(0.0, 0.0, 1.0)
surf4 = TrianglePatch( surf2, 1, 5, "second-tri-patch")
print surf4
surf3.add(surf4)
print surf3

print "Render to VRML"
outfile = open("triangle_patch_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(surf1.vrml_str() + "\n")
outfile.write(surf3.vrml_str() + "\n")
outfile.close()

print "Done."
