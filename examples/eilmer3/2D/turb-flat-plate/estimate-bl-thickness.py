#!/usr/bin/env python
# estimate-bl-thickness.py
#
# Estimates the thickness of the boundary layer by assumming
# boundary-layer-edge as 0.99 of freestream velocity.
#
# Wilson Chan, 11-Jul-2011
#----------------------------------------------------------------

# Filename for file containing boundary layer profile
blProfFileName = "turb_flat_plate-x-368mm.dat"

# Freestream velocity, m/s
U_e = 712.9      # Freestream velocity, m/s
y_wall = 0.16  # y-coordinate of no-slip wall in simulation, m

print "Begin estimate_bl_thickness.py"

print " . Reading boundary layer profile at 0.372 m from file ..", blProfFileName
fi = open(blProfFileName, "r")
fi.readline()  # read away the first line containing notes on file format
#
# Set y-coordinate of the edge of boundary layer to be the same as
# the no-slip wall first.
y_bl_edge = y_wall
#
# Read each line of data and breaking only if we have reached the end of the file,
# or if we have located the edge of the boundary layer.
while True:
    buf = fi.readline().strip()
    if len(buf) == 0: break
    tokens = [float(word) for word in buf.split()]
    # Noting that tokens[1] is the y-coordinate and tokens[5] is the u-velocity,
    # and that turb_flat_plate-x-368mm.dat file is written in an ascending order
    # in the j-th coordinate (so j=0 ... j=1 ... j=n-1), we can locate the boundary
    # layer edge by finding when the u-velocity decreases below 0.99 U_e.
    if tokens[5] <= (0.99 * U_e):
        y_bl_edge = tokens[1]
        break
fi.close()

tol = 1.0e-9
if abs(y_wall - y_bl_edge) < tol:
    print "Boundary layer thickness is less than ", tol, "m."
    print "Warning! The value for boundary layer thickness may be wrong."
else:
    print " . Assuming that edge of boundary layer occurs at 0.99 U_e,"
    print " . where U_e is the freestream velocity of ", U_e , " m/s .." 
    print "boundary_layer_thickness_metres=", abs(y_wall - y_bl_edge)
    print "Done."
