## \file triangle_join_test.py
## \ingroup libgeom2
## \brief Exercise the TrianglePatch Surface class -- joining of patches.
## \author PJ
## \version 05-Mar-2006

from libgeom2 import *

print "\n\n-----------------------------------------------"
print "Begin triangle_join_test..."

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
c0.translate(0.0, 0.0, 1.0)
c1.translate(0.0, 0.0, 1.0)
c2.translate(0.0, 0.0, 1.0)
c3.translate(0.0, 0.0, 1.0)
# surf2 = CoonsPatch(c0, c1, c2, c3, "Shifted-cylinder-surface")
# surf4 = TrianglePatch( surf2, 1, 5, "second-tri-patch")
surf2 = CoonsPatch(c2, c3, c0, c1, "Shifted-cylinder-surface")
surf4 = TrianglePatch( surf2, 5, 1, "second-tri-patch")
print surf4
surf3.add(surf4)
print surf3

print "Render to VRML"
outfile = open("triangle_join_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(surf3.vrml_str() + "\n")
# outfile.write(surf4.vrml_str() + "\n")
outfile.close()

print "Done."
