# test-quad.py
# troublesome example from Rolf's case6 grid.
# PJ, 04-Sep-2012

import sys
import os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
import math
from libprep3 import *
import cfpylib.geom.minimal_geometry as mg

print "Begin using C++ code"

# from running duct6pj calculation, we see these data toward the end.
p0 = Vector3(0.617326, 0.0345716, 1.50002)
p1 = Vector3(0.617389, 0.0345718, 0.500025)
p2 = Vector3(0.618091, 0.0115555, 0.500025)
p3 = Vector3(0.618029, 0.0115554, 1.50002)

centroid = quad_centroid(p0, p1, p2, p3)
n = quad_normal(p0, p1, p2, p3)
t1 = quad_tangent1(p0, p1, p2, p3)
t2 = quad_tangent2(p0, p1, p2, p3)
area = quad_area(p0, p1, p2, p3)

print "centroid=", centroid, "area=", area
print "n=", n, "vabs(n)-1.0=", vabs(n)-1.0
print "t1=", t1, "vabs(t1)-1.0=", vabs(t1)-1.0
print "t2=", t2, "vabs(t2)-1.0=", vabs(t2)-1.0


print "---------------------------------"
print "Begin using pure Python code"

p0=mg.Vector(0.617326, 0.0345716, 1.50002)
p1=mg.Vector(0.617389, 0.0345718, 0.500025)
p2=mg.Vector(0.618091, 0.0115555, 0.500025)
p3=mg.Vector(0.618029, 0.0115554, 1.50002)

centroid, n, t1, t2, area = mg.quad_properties(p0, p1, p2, p3)
print "centroid=", centroid, "area=", area
print "n=", n, "abs(n)-1.0=", abs(n)-1.0
print "t1=", t1, "abs(t1)-1.0=", abs(t1)-1.0
print "t2=", t2, "abs(t2)-1.0=", abs(t2)-1.0

print "Done."
