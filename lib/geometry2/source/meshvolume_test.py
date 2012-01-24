## \file meshvolume_test.py
## \ingroup libgeom2
## \brief Exercise the MeshVolume class.
## \author PJ
## \version 29-Oct-2006

from libgeom2 import *

print "\n\n-----------------------------------------------"
print "Begin meshvolume_test..."

print "\nSimple corner-defined block:"
a = 0.05
p0 = Node(0.0, 0.0, 0.0); p1 = Node(  a, 0.0, 0.0)
p2 = Node(  a,   a, 0.0); p3 = Node(0.0,   a, 0.0)
p4 = Node(0.0, 0.0,   a); p5 = Node(  a, 0.0,   a)
p6 = Node(  a,   a,   a); p7 = Node(0.0,   a,   a)
plist = [p0, p1, p3, p2, p4, p5, p7, p6]
vol1 = MeshVolume(plist, 2, 2, 2, "SIMPLE_BOX")
print vol1
print "vol1.eval(0.25,0.75,0.3)=", vol1.eval(0.25,0.75,0.3)

print "\nBuild volume from mesh that was previously written to a VTK file:"
# Use tim's cylinder as the example
# because it is interesting to mirror-image it.
vol2 = MeshVolume("cyl.0.g")
print vol2
print "vol2.eval(0.25,0.75,0.3)=", vol2.eval(0.25,0.75,0.3)
vol3 = vol2.clone()
vol3.mirror_image(Vector(0.0,0.0,0.0), Vector(0.0,1.0,0.0))

print "\nRender to VRML"
outfile = open("meshvolume_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(vol1.vrml_str() + "\n")
outfile.write(vol2.vrml_str() + "\n")
outfile.write(vol3.vrml_str() + "\n")
outfile.close()

print "Done."
