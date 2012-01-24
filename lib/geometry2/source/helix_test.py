## \file helix_test.py
## \ingroup libgeom2
## \brief Exercise the Helix from the Path C++ classes.
## \author PJ
## \version 24-Apr-2007

from libgeom2 import *
from math import sqrt

print "---------------------------------------------------"
print "Begin helix_test.py ..."

axis0 = Node(0.0, 0.0, 0.0, "a") 
axis1 = Node(1.0, 0.0, 0.0, "b")
pstart = Vector3(0.0, 1.0, 0.0)
pend = Vector3(1.0, 0.0, 1.0)
h1 = Helix(pstart, pend, axis0, axis1, "my_helix")
p = h1.eval(0.5);
print "h1=", h1, "p(0.5)=", p

# make another that crosses the first.
pstart2 = Vector(0.0, 0.0, 1.0)
pend2 = Vector(1.0, 1.0, 0.0)
h2 = Helix(pstart2, pend2, axis0, axis1, "my_helix2", 0.5, 1.0)
p = h2.eval(0.0);
print "h2=", h2, "p(0.0)=", p

print "Render to VRML"
outfile = open("helix_test_py.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(axis0.vrml_str(0.05) + "\n")
outfile.write(axis1.vrml_str(0.05) + "\n")
outfile.write(h1.vrml_str() + "\n")
outfile.write(h2.vrml_str() + "\n")
outfile.close()

print "Done."

