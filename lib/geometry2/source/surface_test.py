## \file surface_test.py
## \ingroup libgeom2
## \brief Exercise the Surface classes.
## \author PJ
## \version 11-Jan-2006

from libgeom2 import *
from math import pi

print "\n\n-----------------------------------------------"
print "Begin surface_test..."

print "\nSimple Coons surface:"
south = Line(Vector3(-1.0, -1.0, 0.0), Vector3(1.0, -0.7, 0.0))
north = Line(Vector3(-1.0, 1.0, 0.0), Vector3(1.0, 1.0, 0.0))
west = Line(Vector3(-1.0, -1.0, 0.0), Vector3(-1.0, 1.0, 0.0))
east = Line(Vector3(1.0, -0.7, 0.0), Vector3(1.0, 1.0, 0.0))
surf1 = AOPatch(south, north, west, east, "AO_SURFACE", 10, 10, 0.25, 0.5, 0.0, 1.0)
print surf1
print "surf1.eval(0.25,0.75)=", surf1.eval(0.25,0.75)

print "\nSimple Cylindrical (Coons) surface:"
c0 = Line(Vector3(1.0, 0.0, 0.0), Vector3(1.0, 0.0, 1.0))
c1 = c0.copy(); c1.translate(-1.0, 1.0, 0.0)
c2 = Arc(Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(0,0,0))
c3 = c2.copy(); c3.translate(0.0, 0.0, 1.0)
surf2 = CoonsPatch(c0, c1, c2, c3, "Cylinder_surface")
print surf2
print "surf2.eval(0.5,0.5)=", surf2.eval(0.5,0.5)
print "surf2.dpds(0.5,0.5)=", surf2.dpds(0.5,0.5), \
      "surf2.dpdr(0.0,0.5)=", surf2.dpdr(0.0,0.5)

print "\nRotate surface 2 about z-axis:"
surf2b = surf2.clone().rotate_about_zaxis(pi/2).rotate_about_zaxis(pi/2)
print "surf2b.eval(0.5,0.5)=", surf2b.eval(0.5,0.5)
print "surf2b.dpds(0.5,0.5)=", surf2b.dpds(0.5,0.5), \
      "surf2b.dpdr(0.0,0.5)=", surf2b.dpdr(0.0,0.5)

print "\nMirror image of the cylindrical surface in the yz-plane."
surf2_mirror = surf2.clone().mirror_image(Vector(0.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0))

print "\nSquare surface patch."
p00 = Node(-1.0, -1.0, 2.0)
p10 = Node(1.0, -1.0, 2.0)
p11 = Node(1.0, 1.0, 2.0)
p01 = Node(-1.0, 1.0, 2.0)
surf3 = CoonsPatch(p00, p10, p11, p01, "Square_surface")
print surf3
print "surf3.eval(0.25,0.75)=", surf3.eval(0.25,0.75)

print "\nPath on the surface (reverse diagonal):"
psurf = PathOnSurface(surf3, LinearFunction(), LinearFunction(-1.0,1.0))
print "psurf=", psurf
n = 5;
dt = 1.0/n
print "    t     position(t)"
for i in range(n+1):
    t = dt * i
    print "pos(", t, ")=", psurf.eval(t) 

print "\nTry a MeshPatch:"
plist = [Vector(0.0, 0.0, 3.0), Vector(0.5, 0.0, 3.5), Vector(1.0, 0.0, 3.0),
         Vector(0.0, 0.8, 3.0), Vector(0.5, 0.8, 3.5), Vector(1.0, 0.8, 3.0)]
surf4 = MeshPatch(plist, 3, 2, "Mesh-surface")
print surf4
print "surf4.eval(0.25,0.75)=", surf4.eval(0.25,0.75)

print "\nTry a Python-powered surface:"
def myfun(r,s): return (2.0*r,3.0*s,0.0)
surf5 = PyFunctionSurface(myfun, "Stretched plane")
print surf5
print "surf5.eval(0.25,0.75)=", surf5.eval(0.25,0.75)

print "\nChannelPatch:"
south = Line(Vector3(-3.0, -1.0, 0.0), Vector3(-2.0, -1.0, 0.0))
north = Arc(Vector3(-3.0, 1.0, 0.0), Vector3(-2.0, 1.0, 0.0), Vector3(-2.5, 1.5, 0.0))
surf6 = ChannelPatch(south, north)
print surf6
print "surf6.eval(0.25,0.75)=", surf6.eval(0.25,0.75)
south.translate(0.0,2.5,0.0)
north.translate(0.0,2.5,0.0)
surf7 = ChannelPatch(south, north, True)
print surf7
print "surf7.eval(0.25,0.75)=", surf7.eval(0.25,0.75)

print "Render to VRML"
outfile = open("surface_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(p00.vrml_str(radius=0.05) + "\n")
outfile.write(p10.vrml_str(radius=0.05) + "\n")
outfile.write(p01.vrml_str(radius=0.05) + "\n")
outfile.write(p11.vrml_str(radius=0.05) + "\n")
outfile.write(surf1.vrml_str() + "\n")
outfile.write(surf2.vrml_str() + "\n")
outfile.write(surf2_mirror.vrml_str() + "\n")
outfile.write(surf3.vrml_str() + "\n")
outfile.write(surf4.vrml_str() + "\n")
outfile.write(surf5.vrml_str() + "\n")
outfile.write(surf6.vrml_str() + "\n")
outfile.write(surf7.vrml_str() + "\n")
outfile.close()

print "Done."
