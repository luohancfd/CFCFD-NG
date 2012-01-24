# reference_temperature_method.py
#
# This script implements the reference temperature method
# described in White (2006), "Viscous fluid flow", 3rd edition,
# pg 517. This method is used to estimate the skin friction 
# for flat plate compressible flow, and is considered to be a
# good approximation in comparison with Van Driest's exact
# calculations.
#
# Wilson Chan, 12 Jan 2010
# PJ refactored it 31-Oct-2010

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from cfpylib.gasdyn import sutherland
from math import sqrt

print "Begin calculating c_f using reference temperature method .."
# Freestream conditions
T_e = 300.0                  # degrees K
p_e = 1.013e3                # static pressure in Pa
u_e = 1390.0                 # u-velocity in m/s
rho_e = p_e / (287.1 * T_e)  # density in kg/m**3
M_e = 4.0                    # u_e / sqrt(1.4 * 287.1 * T_e)
T_w = 300.0                  # Wall temperature in K

outfile = open("cf-ref-temp.data", "w")
outfile.write("# 1.x(m) 2.c_f \n")

x = 0.001 # start just downstream  of the leading edge
x_stop = 1.0
x_step = 0.001
while 1:
    # Reference temperature Eq. 7-42
    T_star = T_e * (0.5 + 0.039 * M_e**2 + 0.5 * T_w / T_e)   # White
    # T_star = T_e * (0.42 + 0.032 * M_e**2 + 0.58 * T_w / T_e)  # Anderson
    c_star = (T_star/T_e)**(-1.0/3.0)   # Eq. 7-40
    # Viscosity by Sutherland's law
    mu_e = sutherland.mu(T_e, "Air")
    Re_xe = rho_e * u_e * x / mu_e     # Eq. 1-1
    c_f = 0.664 * sqrt(c_star) / sqrt(Re_xe)   # Eq. 7-41b
    # c_f = 0.584 / sqrt(Re_xe)
    x += x_step
    outfile.write("%f %f \n" % (x, c_f))
    if x > x_stop: break

outfile.close()
print "Done."
