#!/usr/bin/env python
# estimate_shock_angle.py
# PJ, 11-Jan-2011

import sys, os, gzip
sys.path.append(os.path.expandvars("$HOME/e3bin"))
import math
import numpy
from getopt import getopt
from e3_flow import StructuredGridFlow

#---------------------------------------------------------------
def locate_shock_along_strip(x, y, p):
    """
    Shock location is identified as a pressure rise along a strip of points.

    Input:
    x: sequence of float values, x-coordinate of each point
    y: sequence of float values, y-coordinate
    p: sequence of float values, static pressure in Pa
    Returns:
    x and y coordinates of a point on the shock.

    This function taken from the example 2D/sphere-heat-transfer.
    """
    n = len(x)
    p_max = max(p)
    p_trigger = p[0] + 0.3 * (p_max - p[0])
    x_old = x[0]; y_old = y[0]; p_old = p[0]
    for i in range(1,n):
        x_new = x[i]; y_new = y[i]; p_new = p[i]
        if p_new > p_trigger: break
        x_old = x_new; y_old = y_new; p_old = p_new
    frac = (p_trigger - p_old) / (p_new - p_old)
    x_loc = x_old * (1.0 - frac) + x_new * frac
    y_loc = y_old * (1.0 - frac) + y_new * frac
    return x_loc, y_loc

def locate_shock_front(jobName, nbi, nbj):
    """
    Reads flow blocks and returns the coordinates of the shock front.

    Input:
    jobName: string name used to construct file names
    nbi: number of blocks in the i-index direction
    nbj: number of blocks in the j-index direction

    It is assumed that the shock front will be located by scanning
    along the i-index direction, with j being constant for each search.
    
    This function taken from the example 2D/sphere-heat-transfer
    and has been adapted to the 3D simple_ramp case which has
    the interesting stuff happening in in the x,z-plane.
    """
    blockData = []
    for ib in range(nbi):
        blockData.append([])
        for jb in range(nbj):
            blkindx = ib*nbj + jb
            fileName = 'flow/t9999/%s.flow.b%04d.t9999.gz' % (jobName, blkindx)
            fp = gzip.open(fileName, "r")
            blockData[ib].append(StructuredGridFlow())
            blockData[ib][-1].read(fp)
            fp.close()
    x_shock = []; z_shock = []
    for jb in range(nbj):
        nj = blockData[0][jb].nj
        nk = blockData[0][jb].nk
        for j in range(nj):
            for k in range(nk):
                x = []; z = []; p = [];
                for ib in range(nbi):
                    ni = blockData[ib][jb].ni
                    for i in range(ni):
                        x.append(blockData[ib][jb].data['pos.x'][i,j,k])
                        z.append(blockData[ib][jb].data['pos.z'][i,j,k])
                        p.append(blockData[ib][jb].data['p'][i,j,k])
                xshock, zshock = locate_shock_along_strip(x, z, p)
                x_shock.append(xshock)
                z_shock.append(zshock)
    return x_shock, z_shock

#----------------------------------------------------------------
print "Begin estimate_shock_angle.py"

xs_all, ys_all = locate_shock_front("simple_ramp", nbi=2, nbj=1)
# print "xs_all=", xs_all, "ys=", ys_all
# The shock interacts with the NORTH boundary and 
# the sumulation hasn't reached steady state, so truncate the data.
xs = []; ys = [];
for i in range(len(xs_all)):
    if xs_all[i] < 0.65:
        xs.append(xs_all[i])
        ys.append(ys_all[i])
# print "xs=", xs, "ys=", ys
print "len(xs)=", len(xs), "len(ys)=", len(ys)

# Fit a straight-line to the computed shock points.
m, b = numpy.polyfit(xs, ys, 1) 
print "m=", m, "b=", b
y2 = [m*x+b for x in xs]
shock_angle = math.atan(m)
print "shock_angle_rad=", shock_angle
print "shock_angle_deg=", shock_angle*180/math.pi

# Generate some points on the cone surface.
tan10 = math.tan(10.0*math.pi/180.0)
ycone = [tan10*(x-0.2) for x in xs]

# Average deviation of CFD shock points from fitted line.
d = 0
for i in range(len(xs)):
    d += abs(y2[i] - ys[i])
d /= len(xs)
print "average_deviation_metres=", d

# Optionally do the plot.
if len(sys.argv) > 1 and sys.argv[-1] == "--do-plot":
    import pylab
    pylab.plot(xs, ys, 'o', label='Eilmer3')
    pylab.hold(True)
    pylab.plot(xs, y2, label='fitted line')
    pylab.plot(xs, ycone, label='ramp surface')
    pylab.title('Shock position')
    pylab.xlabel('x, m')
    pylab.ylabel('z, m')
    pylab.legend(loc='upper left')
    pylab.show()

print "Done."
