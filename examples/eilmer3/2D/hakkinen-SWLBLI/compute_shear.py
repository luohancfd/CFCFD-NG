#! /usr/bin/env python
# compute_shear.py
#
# Pick up the simulation data at the last simulated time
# and compute an estimate of the shear-stress coefficient.
#
# PJ, 08-May-2013

import sys, os
from math import sqrt
from e3_flow import read_all_blocks
job = "swlbli"
nb = 20
pick_list = [2, 4, 6, 8, 10, 12, 14, 16, 18] # blocks against plate
rho_inf = 0.1315 # kg/m**3
u_inf =  514.0 # m/s
T_inf = 164.4 # K
from cfpylib.gasdyn import sutherland
mu_inf = sutherland.mu(T_inf, 'Air')

print "Determine the latest time."
fp = open(job+".times", "r"); lines = fp.readlines(); fp.close()
tindx = int(lines[-1].strip().split()[0]) # first number of the last line
print "tindx=", tindx
print "Begin: Pick up data."
grid, flow, dim = read_all_blocks(job, nb, tindx, zipFiles=True)
print "Compute shear stress for cell-centres along plate surface"
outfile = open("shear.data", "w")
outfile.write("# x(m)  tau_w(Pa)  Cf   y_plus\n")
for ib in pick_list:
    j = 0 # plate is along the South boundary  
    k = 0 # of a 2D grid
    print "# start of block"
    for i in range(flow[ib].ni):
        # Cell closest to surface
        x = flow[ib].data['pos.x'][i,j,k]; 
        y = flow[ib].data['pos.y'][i,j,k]
        rho = flow[ib].data['rho'][i,j,k]; 
        u1 = flow[ib].data['vel.x'][i,j,k]
        mu = flow[ib].data['mu'][i,j,k]
        dudy = (u1 - 0.0) / y # Assuming that the wall is straight down at y=0
        tau_w = mu * dudy    # wall shear stress
        Cf = tau_w / (0.5*rho_inf*u_inf*u_inf)
        u_tau = sqrt(abs(tau_w) / rho) # friction velocity
        y_plus = u_tau * y * rho / mu
        Rex = rho_inf * u_inf * x / mu_inf
        Cf_blasius = 0.664 / sqrt(Rex) 
        outfile.write("%f %f %f %f %f\n" % (x, tau_w, Cf, Cf_blasius, y_plus))
        print "x=", x, "tau_w=", tau_w, "Cf=", Cf, "y_plus=", y_plus
outfile.close()
print "Done"
