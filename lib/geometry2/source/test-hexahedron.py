# test-hexahedron.py
# troublesome example from  Rolf and Rowan's duct-inject case
# PJ, 10-Sep-2012

import sys
import os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
import math
from libprep3 import *
import cfpylib.geom.minimal_geometry as mg

if 1:
    # data copied from LOGFILE-1
    # hexahedron_properties(): significant negative volume: -7.56648776166e-07
    fy = -1.0 # for flipping cell over
    p0=Vector3(1.2684011193, -0.000915290907604*fy, 0.497645306831) 
    p1=Vector3(1.26901312421, -0.000558656181196*fy, 0.497989150727) 
    p2=Vector3(1.27283724882, -0.00617661471085*fy, 0.452538471634) 
    p3=Vector3(1.2727171364, -0.00562783086504*fy, 0.452361419373)
    p0, p1 = p1, p0  # I think that the original data is twisted.
    p4=Vector3(1.21921640281, -0.00702801321524*fy, 0.516907489331) 
    p5=Vector3(1.2192610028, -0.00767192831942*fy, 0.517088081703) 
    p6=Vector3(1.22574584479, -0.0220765393669*fy, 0.473315349638) 
    p7=Vector3(1.226594613, -0.0212290198649*fy, 0.473503835374)
else:
    p0 = Vector3(0.0, 0.0, 0.0)
    p1 = Vector3(1.0, 0.0, 0.0)
    p2 = Vector3(1.0, 1.0, 0.0)
    p3 = Vector3(0.0, 1.0, 0.0)
    p4 = Vector3(0.0, 0.0, 0.01)
    p5 = Vector3(1.0, 0.0, 0.01)
    p6 = Vector3(1.0, 1.0, 0.01)
    p7 = Vector3(0.0, 1.0, 0.01)

print "Begin using C++ hexahedron volume code"
vol = hexahedron_volume(p0, p1, p2, p3, p4, p5, p6, p7)
cent = hexahedron_centroid(p0, p1, p2, p3, p4, p5, p6, p7)
print "vol=", vol, "centroid=", cent
print "-----------------------------"

print "Begin using C++ hexahedron2 volume code"
vol = hexahedron2_volume(p0, p1, p2, p3, p4, p5, p6, p7)
cent = hexahedron2_centroid(p0, p1, p2, p3, p4, p5, p6, p7)
print "vol=", vol, "centroid=", cent
print "-----------------------------"

print "Begin using C++ hex_cell volume code"
vol = hex_cell_volume(p0, p1, p2, p3, p4, p5, p6, p7)
cent = hex_cell_centroid(p0, p1, p2, p3, p4, p5, p6, p7)
print "vol=", vol, "centroid=", cent

vol1 = SimpleBoxVolume(p0, p1, p2, p3, p4, p5, p6, p7, "SIMPLE_BOX")
print "\nRender to VRML"
outfile = open("test-hexahedron.wrl", "w")
outfile.write("#VRML V2.0 utf8\n")
outfile.write(vol1.vrml_str() + "\n")
outfile.close()
print "Done."
