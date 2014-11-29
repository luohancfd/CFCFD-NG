# compute_viscous_data_simple.py
# inspired by PJ's compute_y_plus.py script
# Updated 27 October 2008 - Wilson Chan
# Refined 01 November 2010 - PJ
# Altered for 3D Flat Plate May 2014 - Samuel Stennett

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from cfpylib.gasdyn import sutherland

# Freestream and test conditions (user input)
wall_y = 0.16      # Location of wall, m
T_w = 316.2        # Temperature at wall, K
rho_inf = 0.177    # free stream density, kg/m^3
u_inf = 712.89     # free stream velocity, m/s
mu = sutherland.mu(T_w, "Air")  # fluid viscosity at wall

print "Begin compute_viscous_data_master.py..."

#SLICE 0

print "Begin compute_viscous_data_simple.py for slice 0..."

infile = open("./turb_flat_plate-y-wall0.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data0.dat", "w")
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

#SLICE 1

print "Begin compute_viscous_data_simple.py for slice 1..."

infile = open("./turb_flat_plate-y-wall1.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data1.dat", "w")
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

#SLICE 2

print "Begin compute_viscous_data_simple.py for slice 2..."

infile = open("./turb_flat_plate-y-wall2.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data2.dat", "w")
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

#SLICE 3

print "Begin compute_viscous_data_simple.py for slice 3..."

infile = open("./turb_flat_plate-y-wall3.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data3.dat", "w")
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

#SLICE 4

print "Begin compute_viscous_data_simple.py for slice 4..."

infile = open("./turb_flat_plate-y-wall4.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data4.dat", "w")
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

#SLICE 5

print "Begin compute_viscous_data_simple.py for slice 5..."

infile = open("./turb_flat_plate-y-wall5.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data5.dat", "w")
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

#SLICE 6

print "Begin compute_viscous_data_simple.py for slice 6..."

infile = open("./turb_flat_plate-y-wall6.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data6.dat", "w")
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

#SLICE 7

print "Begin compute_viscous_data_simple.py for slice 7..."

infile = open("./turb_flat_plate-y-wall7.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data7.dat", "w")
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

#SLICE 8

print "Begin compute_viscous_data_simple.py for slice 8..."

infile = open("./turb_flat_plate-y-wall8.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data8.dat", "w")
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

#SLICE 9

print "Begin compute_viscous_data_simple.py for slice 9..."

infile = open("./turb_flat_plate-y-wall9.dat", "r")
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("viscous_data9.dat", "w")
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

print "compute_viscous_data_master.py complete!"
