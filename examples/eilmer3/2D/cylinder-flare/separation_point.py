#! /usr/bin/env python
# separation_point.py
#
# Pick up the simulation data at all time frames.
# Search for the zero-crossing of ux to identify the separation point
# on the cylinder surface.
#
# PJ, 07-June-2013
print "Begin..."
import sys, os
from e3_flow import read_all_blocks
#
nb = 24
pick_list = [0, 2, 4, 6, 8, 10] # blocks against cylinder only
job = "cyl-flare"
fp = open(job+".times", "r"); lines = fp.readlines(); fp.close()
times = []; xzero = []
for item in lines:
    items = item.strip().split()
    if items[0] == '#': continue
    tindx = int(items[0])
    if tindx == 0: continue
    t = float(items[1])
    print "Begin: Pick up data for tindx=", tindx, "t=", t
    grid, flow, dim = read_all_blocks(job, nb, tindx, zipFiles=True)
    x = []; y = []; ux = []
    for ib in pick_list:
        j = 0 # surface is along the South boundary  
        k = 0 # of a 2D grid
        for i in range(flow[ib].ni):
            # Cell closest to surface
            x.append(flow[ib].data['pos.x'][i,j,k]) 
            ux.append(flow[ib].data['vel.x'][i,j,k])
    # Find the zero-crossing interval, 
    # assuming that we start with positive velocity.
    # For no zero-crossing we run to the end.
    i = 0
    while ux[i] >= 0.0 and i < len(ux)-1: i += 1
    # Linearly interpolate the zero-crossing point.
    frac = ux[i-1]/(ux[i-1]-ux[i])
    xzero.append((1.0-frac)*x[i-1] + frac*x[i])
    times.append(t)
    print "t=", t, "xzero=", xzero[-1]

outfile = open("separation-location.data", "w")
outfile.write("# t(s) x(m)\n")
for i in range(len(xzero)):
    outfile.write("%f %f\n" % (times[i], xzero[i]))
outfile.close()

outfile = open("separation-velocity.data", "w")
outfile.write("# t(s) -dx/dt(m/s)\n")
for i in range(1,len(xzero)):
    outfile.write("%f %f\n" % (times[i], -(xzero[i]-xzero[i-1])/(times[i]-times[i-1])))
outfile.close()

print "Fit an asymptotic function to the location data."
import numpy
x = numpy.array(xzero) * 1000.0 # to get units of mm
t = numpy.array(times) * 1000.0 # to get units of ms
def f(t, xf, dx, tau):
    return xf + dx * numpy.exp(-t/tau)
from scipy.optimize import curve_fit
popt, pcov = curve_fit(f, t, x, [60.0, 30.0, 0.8])
print "Fitted parameters:"
print "xf=", popt[0], "mm"
print "dx=", popt[1], "mm"
print "tau=", popt[2], "ms"
print "pcov=", pcov
print "Done"
