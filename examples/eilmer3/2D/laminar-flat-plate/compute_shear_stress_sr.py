# compute_shear_stress_sr.py
# Blame: PJ and Wilson Chan.
# 
import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from math import sqrt
from e3_grid import *
from e3_block import *
from e3_flow import *
from cfpylib.gasdyn import sutherland

print "Begin shear_stress.py..."
T_w = 300.0  # Wall temperature in degrees K
# Freestream values of density and x-velocity
# to nondimensionalise the shear-stress data.
rho_inf = 0.011764 # kg/m^3
u_inf = 1390.0     # m/s

print "Pick up data..."
grid, flow, dimensions = read_all_blocks(rootName="lam_flat_plate_sr", 
                                         nblock=16, tindx=9999, zipFiles=1)

print "Compute viscous parameters for cell-centres along plate surface:"
outfile = open("cf-eilmer3-sr.data", "w")
outfile.write("# 1.x(m) 2.du(m/s) 3.dy(m) 4.tau_w(N/m^2) 5.c_f 6.y_plus \n")

wall_block_index_list = [3,7,11,15] # along the NORTH wall.
for count in range(len(wall_block_index_list)):
    ib = wall_block_index_list[count]
    j = flow[ib].nj - 1 # first cell from wall
    for i in range(flow[ib].ni):
        # Work along the plate for the current block.
        x = flow[ib].data['pos.x'][i][j]; y = flow[ib].data['pos.y'][i][j]
        # Get density & u-velocity of cell
        rho = flow[ib].data['rho'][i][j]
        u1 = flow[ib].data['vel.x'][i][j]
        # Compute distance from centre to interface (2D only for now)
        c = Vector(x[0], y[0])  # cell centre
        # north-west vertex
        temp_x = grid[ib].x[i,j+1]; temp_y = grid[ib].y[i,j+1] 
        s0 = Vector(temp_x[0], temp_y[0])
        # north-east vertex
        temp_x = grid[ib].x[i+1,j+1]; temp_y = grid[ib].y[i+1,j+1]
        s1 = Vector(temp_x[0], temp_y[0])
        tangent = unit(s1 - s0)
        s0c = c - s0
        # d1 is distance from centre to interface
        d1 = vabs(s0c - dot(s0c,tangent)*tangent)
        # Compute viscous parameters
        dudy = (u1 - 0.0) / d1   # du/dy for no-slip wall
        mu = sutherland.mu(T_w, "Air")
        # Wall shear stress .. tau_w = mu_w * du/dy
        tau_w = mu * dudy
        # Friction velocity .. u_tau = (tau_w / rho_w)**0.5
        u_tau = sqrt(tau_w / rho)
        # Non-dimensionalised y .. y_plus = u_tau * dy * rho_w / mu_w
        y_plus = u_tau * d1 * rho / mu
        # Skin friction coefficient .. c_f = tau_w / (0.5 * rho * u**2)
        c_f = tau_w / (0.5 * rho_inf * u_inf**2)
        outfile.write("%f %f %f %f %f %f \n" % (x, u1, d1, tau_w, c_f, y_plus))

outfile.close()
print "Done."

