#! /usr/bin/env python
# surface_properties.py
#
# Pick up the simulation data at the last simulated time
# compute an estimate of the shear-stress coefficient and
# output both shear and pressure along the cubic surface.
#
# PJ, 11-Aug-2013

import sys, os
job = "cubic-ramp"
print "Determine the latest time."
fp = open(job+".times", "r"); lines = fp.readlines(); fp.close()
tindx = int(lines[-1].strip().split()[0]) # first number of the last line
print "tindx=", tindx

print "Begin: Pick up data for tindx=", tindx
from libprep3 import Vector, cross, dot, vabs
from e3_flow import read_all_blocks
from math import sqrt
#
nb = 22
pick_list = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18] # blocks against cubic surface
rho_inf = 5.521e-3 # kg/m**3
p_inf = 66.43 # Pa
u_inf =  1589.8 # m/s
T_inf = 41.92 # K
T_wall = 296.0 # K
from cfpylib.gasdyn import sutherland
mu_inf = sutherland.mu(T_inf, 'Air')
mm = 0.001  # metres
#
grid, flow, dim = read_all_blocks(job, nb, tindx, zipFiles=True)
print "Compute shear stress for cell-centres along the surface"
outfile = open("surface.data", "w")
outfile.write("# x(m) tau_w(Pa) Cf Cf_blasius y_plus p(Pa) Cp q(W/m**2) Ch\n")
for ib in pick_list:
    j = 0 # surface is along the South boundary  
    k = 0 # of a 2D grid
    print "# start of block"
    for i in range(flow[ib].ni):
        # Cell closest to surface
        x = flow[ib].data['pos.x'][i,j,k] 
        y = flow[ib].data['pos.y'][i,j,k]
        ctr = Vector(x, y)
        # Get vertices on surface, for this cell.
        x = grid[ib].x[i,j,k] 
        y = grid[ib].y[i,j,k]
        vtx0 = Vector(x, y)
        x = grid[ib].x[i+1,j,k] 
        y = grid[ib].y[i+1,j,k]
        vtx1 = Vector(x, y)
        t1 = (vtx1-vtx0)
        t1.norm() # tangent vector for surface
        midpoint = 0.5*(vtx0+vtx1) # on surface
        normal = cross(Vector(0,0,1),t1)
        normal.norm()
        # Surface to cell-centre distance.
        dy = dot(normal, ctr-midpoint)
        # Cell-centre flow data.
        rho = flow[ib].data['rho'][i,j,k] 
        ux = flow[ib].data['vel.x'][i,j,k]
        uy = flow[ib].data['vel.y'][i,j,k]
        v = Vector(ux, uy)
        vt = dot(v,t1) # velocity component tangent to surface
        mu = flow[ib].data['mu'][i,j,k]
        kgas = flow[ib].data['k[0]'][i,j,k]
        p = flow[ib].data['p'][i,j,k]
        Cp = (p-p_inf)/(0.5*rho_inf*u_inf*u_inf)
        T = flow[ib].data['T[0]'][i,j,k]
        # Shear stress
        dudy = (vt - 0.0) / dy # no-slip wall
        tau_w = mu * dudy    # wall shear stress
        Cf = tau_w / (0.5*rho_inf*u_inf*u_inf)
        u_tau = sqrt(abs(tau_w) / rho) # friction velocity
        y_plus = u_tau * dy * rho / mu
        Rex = rho_inf * u_inf * midpoint.x / mu_inf
        Cf_blasius = 0.664 / sqrt(Rex) 
        # Heat flux
        dTdy = (T - T_wall) / dy # conductive heat flux at the wall
        q = kgas * dTdy
        Ch = q / (0.5*rho_inf*u_inf*u_inf*u_inf)
        #
        outfile.write("%f %f %f %f %f %f %f %f %f\n" % 
                      (midpoint.x, tau_w, Cf, Cf_blasius,
                       y_plus, p, Cp, q, Ch))
        print "x=", midpoint.x, "tau_w=", tau_w, "Cf=", Cf, "y_plus=", y_plus, \
            "p=", p, "Cp=", Cp, "q=", q, "Ch=", Ch
outfile.close()
print "Done"
