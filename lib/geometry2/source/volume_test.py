## \file volume_test.py
## \ingroup libgeom2
## \brief Exercise the ParametricVolume classes.
## \author PJ
## \version 13-Jan-2006

from libgeom2 import *
from math import pi

print "\n\n-----------------------------------------------"
print "Begin volume_test..."

print "\nSimple corner-defined block:"
a = 1.0
p0 = Node(0.0, 0.0, 0.0); p1 = Node(  a, 0.0, 0.0)
p2 = Node(  a,   a, 0.0); p3 = Node(0.0,   a, 0.0)
p4 = Node(0.0, 0.0,   a); p5 = Node(  a, 0.0,   a)
p6 = Node(  a,   a,   a); p7 = Node(0.0,   a,   a)
vol1 = SimpleBoxVolume(p0, p1, p2, p3, p4, p5, p6, p7, "SIMPLE_BOX",
                       0.0, 0.5, 0.0, 0.5, 0.5, 1.0)
print vol1
print "vol1.eval(0.25,0.75,0.3)=", vol1.eval(0.25,0.75,0.3)

print "\nWire-frame definition for part cylinder block:"
r1 = 0.5; r2 = 0.9;
zero = Vector(0.0, 0.0, 0.0)
p0 = Vector(0.0, r1, 0.0); p1 = Vector(r1, 0.0, 0.0)
p3 = Vector(0.0, r2, 0.0); p2 = Vector(r2, 0.0, 0.0)
c01 = Arc(p0,p1,zero)
c32 = Arc(p3,p2,zero)
c03 = Line(p0,p3)
c12 = Line(p1,p2)
up = Vector(0.0, 0.0, 0.6)
upright = Line(zero, up)
c04 = upright.copy().translate(p0)
print "upright=", upright, "c04=", c04
c15 = upright.copy(); c15.translate(p1)
c26 = upright.copy(); c26.translate(p2)
c37 = upright.copy(); c37.translate(p3)
c45 = c01.copy(); c45.translate(up);
c56 = c12.copy(); c56.translate(up);
c76 = c32.copy(); c76.translate(up);
c47 = c03.copy(); c47.translate(up);
vol2 = WireFrameVolume(c01, c12, c32, c03, c45, c56, c76, c47, c04, c15, c26, c37)
print vol2
print "vol2.eval(0.25,0.75,0.3)=", vol2.eval(0.25,0.75,0.3)

print "\nRotate bolck about the z-axis:"
vol2a = vol2.clone()
vol2a.rotate_about_zaxis(pi)
print "vol2a.eval(0.25,0.75,0.3)=", vol2a.eval(0.25,0.75,0.3)

print "\nCopy and translate a block."
vol3 = vol2.copy()
vol3.translate(1.0, 0.0, 0.0)
print vol3
print "vol3.eval(0.25,0.75,0.3)=", vol3.eval(0.25,0.75,0.3)

print "\nCopy block and mirror in xz-plane."
vol3a = vol2.copy()
vol3a.mirror_image(Vector(0.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0))
print vol3a
print "vol3a.eval(0.25,0.75,0.3)=", vol3a.eval(0.25,0.75,0.3)

print "\nExtruded surface."
base = CoonsPatch( c01, c32, c03, c12 )
upright.translate(Vector(0.0, 0.0, 1.5))
print "upright=", upright
vol4 = WireFrameVolume(base, upright)
print "vol4.eval(0.25,0.75,0.3)=", vol4.eval(0.25,0.75,0.3)

print "\nSurfaceThruVolume."
fr = BilinearFunction(0.0,1.0,0.0,1.0)
fs = BilinearFunction(0.0,0.0,1.0,1.0)
ft = BilinearFunction(0.0,0.0,0.0,0.0)
surf_bottom = SurfaceThruVolume(vol4, fr, fs, ft)
print "surf_bottom=", surf_bottom
print "surf_bottom.eval(0.25,0.75)=", surf_bottom.eval(0.25,0.75)

print "\nPyFunctionVolume:"
def myfun(r,s,t): return(-2.0*r,-3.0*s,-3.0*t)
vol5 = PyFunctionVolume(myfun, "Box-from-myfun")
print "vol5=", vol5
print "vol5.eval(0.5,0.5,0.5)=", vol5.eval(0.5,0.5,0.5)

print "\nRender to VRML"
outfile = open("volume_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(vol1.vrml_str() + "\n")
outfile.write(vol2.vrml_str() + "\n")
outfile.write(vol3.vrml_str() + "\n")
outfile.write(vol3a.vrml_str() + "\n")
outfile.write(vol4.vrml_str() + "\n")
outfile.write(vol5.vrml_str() + "\n")
outfile.close()

print "Done."
