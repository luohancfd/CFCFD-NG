#! /usr/bin/env python
# compute_y_plus.py
#
# Pick up the simulation data at t=5ms for the Coles flat-plate simulation,
# locate the appropriate x-location for writing boundary-layer profiles,
# and compute y+ for the cells nearest the plate.
#
# PJ, 26-Nov-2007

import sys, os
sys.path.append(os.path.expandvars("$HOME/cfd_bin"))
from cfpylib.flow.blockflow2d import BlockFlow2D
from blockgrid2d import BlockGrid2D
from libgeom2 import *

#---------------------------------------------------------------------
print "Begin: Pick up data..."

nblock = 4
t_target = 5.0e-3

print "Read grid data..."
gridfile = open("coles.g", "r")
flow = []; grid = []
for i in range(nblock):
    grid.append(BlockGrid2D())
    flow.append(BlockFlow2D(k_omega=1))
    grid[i].read(gridfile)
    print "Have read grid of", grid[i].ni, "by", grid[i].nj, "vertices."
gridfile.close()

print "Read flow data..."
flowfile = open("coles.s", "r")
flow[0].t = 0.0
while flow[0].t < t_target:
    for i in range(nblock): flow[i].read(flowfile)
flowfile.close()
print "Have read data for time t=", flow[0].t

#---------------------------------------------------------------------
print "Locate index for profile location."
x_location = 0.546
ib_closest = 0
i_closest = 0
dx_closest = 10.0 # something unreasonably large
for ib in range(nblock):
    for i in range(flow[ib].ni):
        dx = abs(flow[ib].x[i][0] - x_location)
        if dx < dx_closest:
            ib_closest = ib
            i_closest = i
            dx_closest = dx
print "Indices of closest wall cell: ib=", ib_closest, \
      "i=", i_closest, \
      "dx=", dx_closest

#--------------------------------------------------------------------
print "Compute y+ for cell-centres along plate surface:"

def viscosity(T):
    "Compute viscosity for air."
    S  = 110.4   # Sutherland constant
    mu_ref = 17.89e-6  # Pa.s
    T_ref = 273.1      # degrees K
    T_T0 = T / T_ref
    mu = mu_ref * T_T0 * sqrt(T_T0)*(T_ref + S)/(T + S)
    return mu

T_w = 287.9  # degrees K
from math import sqrt
outfile = open("y_plus.dat", "w")
outfile.write("# x(metres) tau_w(Pa) y_plus\n")

for ib in range(nblock):
    j = flow[ib].nj - 1 # plate is along the NORTH boundary
    print "# start of block"
    outfile.write("# start of block\n")
    for i in range(flow[ib].ni):
        # Cell closest to surface
        x = flow[ib].x[i][j]; y = flow[ib].y[i][j]
        rho = flow[ib].rho[i][j]
        c = Vector(x, y)  # cell centre
        s0 = Vector(grid[ib].x[i][j+1], grid[ib].y[i][j+1])  # north-west vertex
        s1 = Vector(grid[ib].x[i+1][j+1], grid[ib].y[i+1][j+1]) # north-east vertex
        tangent = unit(s1 - s0) 
        s0c = c - s0
        d1 = vabs(s0c - dot(s0c,tangent)*tangent) # distance from centre to interface
        u1 = flow[ib].vx[i][j]
        dudy = (u1 - 0.0) / d1
        mu = viscosity(T_w)  # viscosity of the gas at the wall
        tau_w = mu * dudy    # wall shear stress
        u_tau = sqrt(tau_w / rho) # friction velocity
        y_plus = u_tau * d1 * rho / mu
        outfile.write("%f %f %f\n" % (x, tau_w, y_plus))
        print "x=", x, "tau_w=", tau_w, "y_plus=", y_plus

outfile.close()
print "Done"
