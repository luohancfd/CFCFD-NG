#! /usr/bin/env python
"""
locate_shock.py -- Locate the shock by its pressure jump.

PJ, 12-Apr-2012
"""
print "Begin..."
import sys, os, gzip
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_flow import StructuredGridFlow

# Block 1 contains the shock and the fully-expanded driver gas.
fileName = 'flow/t9999/cst.flow.b0001.t9999.gz'
fp = gzip.open(fileName, "r")
blockData = StructuredGridFlow()
blockData.read(fp)
fp.close()

# We expect the shock to have progressed some way along the i-index.
# Start the search from the right and move left.
k = 0; j = 0; i = blockData.ni -1
p_trigger = 2.0e6  # Pa
x_old = blockData.data['pos.x'][i,j,k]
p_old = blockData.data['p'][i,j,k]
while i >= 0:
    i -= 1
    x = blockData.data['pos.x'][i,j,k]
    p = blockData.data['p'][i,j,k]
    if p > p_trigger: break
    x_old, p_old = x, p
frac = (p_trigger - p_old) / (p - p_old)
x_loc = x_old * (1.0 - frac) + x * frac
t_final = 100.0e-6 # seconds
print "shock at x=", x_loc, "m, speed=", (x_loc - 0.5)/t_final, "m/s"

# Also compute average gas speed of the expanded driver gas
# over a representative region.
u_sum = 0; n = 0;
for i in range(blockData.ni):
    x = blockData.data['pos.x'][i,j,k]
    u = blockData.data['vel.x'][i,j,k]
    if x >= 0.7 and x <= 0.8:
        u_sum += u; n += 1
u_sum /= n
print "average u_g=", u_sum, "m/s"
print "Done."
