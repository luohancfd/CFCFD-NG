# compute_viscous_data_simple.py
# inspired by PJ's compute_y_plus.py script
# Updated 27 October 2008 - Wilson Chan
# Refined 01 November 2010 - PJ

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from cfpylib.gasdyn import sutherland

# Freestream and test conditions (user input)
wall_y = 0.0      # Location of wall, m
T_w = 316.2        # Temperature at wall, K
rho_inf = 0.457    # free stream density, kg/m^3
u_inf = 707.0     # free stream velocity, m/s
mu = sutherland.mu(T_w, "Air")  # fluid viscosity at wall

print "Begin compute_viscous_data_simple.py"

infile = open("./tc2update-y-wall.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data.dat", "w")
outfile.write("# x(m) du(m/s) dy(m) tau_w(N/m^2) c_f y_plus \n")

# Process every data line.
while True:
    line = infile.readline().strip()
    if len(line) == 0: break
    data_array = [float(word) for word in line.split()]
    x = data_array[0]
    dy = abs(data_array[2] - wall_y)
    du = abs(data_array[5] - 0.0)
    rho = data_array[4]
    mu = data_array[10]
    tau_w = mu * du/dy   # Wall shear stress  
    cf = tau_w / (0.5 * rho_inf * u_inf**2)   # Skin friction coefficient
    u_tau = (tau_w / rho)**0.5   # Friction velocity
    yplus = u_tau * dy * rho / mu
    outfile.write("%f %f %e %f %f %f \n" % (x, du, dy, tau_w, cf, yplus))
    # Output for Eilmer3 automated tests - prints skin friction coefficient at x = 0.366559 m.
    if abs(x - 0.366559) < 1.0e-7:
        print "skin_friction_coefficient_at_367mm=", cf
infile.close()
outfile.close()

print "Done."
