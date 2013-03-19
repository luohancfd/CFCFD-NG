#! /usr/bin/env python
"""
estimate_flow_properties.py

Locate the shock by its pressure jump and estimate shocked
and expanded flow conditions.

PJ, 19-Mar-2013
Adapted from the locate_shock.py script in classic_shock_tube case.
"""
print "Begin..."
import sys, os, gzip
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_flow import StructuredGridFlow

# Block 1 contains the shock and the fully-expanded driver gas.
fileName = 'flow/t9999/sod.flow.b0001.t9999.gz'
fp = gzip.open(fileName, "r")
blockData = StructuredGridFlow()
blockData.read(fp)
fp.close()

# We expect the shock to have progressed some way along the i-index.
# Start the search from the right and move left.
k = 0; j = 0; i = blockData.ni -1
p_trigger = 20000.0  # Pa
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
t_final = 0.6063e-3 # seconds
print "shock at x=", x_loc, "m, speed=", (x_loc - 0.5)/t_final, "m/s"

# Compute average gas speed of the expanded driver gas
# over a representative region.
u_sum = 0; p_sum = 0; e_sum = 0; T_sum = 0; n = 0
for i in range(blockData.ni):
    x = blockData.data['pos.x'][i,j,k]
    if x >= 0.500 and x <= 0.631:
        u_sum += blockData.data['vel.x'][i,j,k]
        p_sum += blockData.data['p'][i,j,k]
        e_sum += blockData.data['e[0]'][i,j,k]
        T_sum += blockData.data['T[0]'][i,j,k]
        n += 1
u_sum /= n; p_sum /= n; e_sum /= n; T_sum /= n
print "expanded-driver-gas averages for n=", n
print "    edg u_g=", u_sum, "m/s"
print "    edg p=", p_sum, "Pa"
print "    edg e=", e_sum, "J/kg"
print "    edg T=", T_sum, "K"

# Compute average gas speed of the shock-compressed driven gas
# over a representative region.
u_sum = 0; p_sum = 0; e_sum = 0; T_sum = 0; n = 0
for i in range(blockData.ni):
    x = blockData.data['pos.x'][i,j,k]
    if x >= 0.723 and x <= 0.810:
        u_sum += blockData.data['vel.x'][i,j,k]
        p_sum += blockData.data['p'][i,j,k]
        e_sum += blockData.data['e[0]'][i,j,k]
        T_sum += blockData.data['T[0]'][i,j,k]
        n += 1
u_sum /= n; p_sum /= n; e_sum /= n; T_sum /= n
print "shocked-driven-gas averages for n=", n
print "    sdg u_g=", u_sum, "m/s"
print "    sdg p=", p_sum, "Pa"
print "    sdg e=", e_sum, "J/kg"
print "    sdg T=", T_sum, "K"
print "Done."
