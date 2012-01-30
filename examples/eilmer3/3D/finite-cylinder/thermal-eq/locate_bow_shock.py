#! /usr/bin/env python
# \file locate_bow_shock.py
# PJ, 08-Nov-2009, updated for Eilmer3

import sys, os, gzip
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_flow import StructuredGridFlow

print "Locate a bow shock by its pressure jump."

# Block 0 contains the stagnation point.
fileName = 'flow/t9999/cyl.flow.b0000.t9999.gz'
fp = gzip.open(fileName, "r")
blockData = StructuredGridFlow()
blockData.read(fp)
fp.close()

# Since this is a 3D simulation, the shock is not expected
# to be flat in the k-direction (along the cylinder axis).
# Sample the shock layer in a few places near the stagnation line.
x_sum = 0.0
n_sample = 6
for k in range(n_sample):
    j = 0
    p_trigger = 10000.0  # Pa
    x_old = blockData.data['pos.x'][0,j,k]
    p_old = blockData.data['p'][0,j,k]
    for i in range(blockData.ni):
        x = blockData.data['pos.x'][i,j,k]
        p = blockData.data['p'][i,j,k]
        if p > p_trigger: break
        x_old = x
        p_old = p
    frac = (p_trigger - p_old) / (p - p_old)
    x_loc = x_old * (1.0 - frac) + x * frac
    print "shock at x=", x_loc, \
        "y=", blockData.data['pos.y'][0,j,k], \
        "z=", blockData.data['pos.z'][0,j,k]
    x_sum += x_loc
x_average = x_sum / n_sample
print "Average x-location=", x_average

print "Done."
