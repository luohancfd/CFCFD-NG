# compute_cp.py
# modified from WC,PJ's compute_viscous_data_simple.py script
# Modified 4 December 2014 by Samuel Stennett

import sys, os, math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from cfpylib.gasdyn import sutherland

# Freestream and test conditions (user input)
P_inf = 7100.0			#Freestream Pressure, Pa
T_inf = 70.3        		#Freestream Temperature, K
rho_inf = P_inf/(287.1*T_inf)   #Freestream Density, kg/m^3
vel_inf = 4.0*math.sqrt(1.4*287.0*T_inf)

#jetD = 0.00412 # OR 0.00476, need to check results (CFD vs actual diameter)
jetD = 0.00476
jetdist = 0.07826

print "Begin compute_cp_normalised.py"

infile = open("./inject-y-slice-UPSTREAM-31.dat", "r") #NEED TO ENSURE THIS IS NAMED CORRECTLY
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("cp-centreline-UPSTREAM-31.dat", "w")
outfile.write("# dist-to-jet(mm) x/D() c_p \n")

# Process every data line.
while True:
    line = infile.readline().strip()
    if len(line) == 0: break
    data_array = [float(word) for word in line.split()]
    surfacepressure = data_array[8]
    currentpos = data_array[0]
    cpval = (surfacepressure-P_inf)/(0.5*rho_inf*vel_inf**2)
    disttojet = currentpos - jetdist
    XonD = disttojet/jetD

    outfile.write("%f %f %f \n" % (disttojet, XonD, cpval))

infile.close()
outfile.close()

print "Done... Calculating downstream data..."

infile = open("./inject-y-slice-DOWNSTREAM-31.dat", "r") #NEED TO ENSURE THIS IS NAMED CORRECTLY
line = infile.readline()  # read away the first (comment) line containing notes on format
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("cp-centreline-DOWNSTREAM-31.dat", "w")
outfile.write("# dist-to-jet(mm) x/D() c_p \n")

# Process every data line.
while True:
    line = infile.readline().strip()
    if len(line) == 0: break
    data_array = [float(word) for word in line.split()]
    surfacepressure = data_array[8]
    currentpos = data_array[0]
    cpval = (surfacepressure-P_inf)/(0.5*rho_inf*vel_inf**2)
    disttojet = currentpos - jetdist
    XonD = disttojet/jetD

    outfile.write("%f %f %f \n" % (disttojet, XonD, cpval))

infile.close()
outfile.close()

print "Combining solutions into one .dat file..."

infile1 = open("./cp-centreline-UPSTREAM-31.dat", "r")
infile2 = open("./cp-centreline-DOWNSTREAM-31.dat", "r")
line1 = infile1.readline()
line2 = infile2.readline()

outfile = open("cp-centreline-complete-31.dat", "w")
outfile.write("# dist-to-jet(mm) x/D() c_p \n")

while True:
    line1 = infile1.readline()
    if len(line1) == 0: 
	break
    outfile.write(line1)

outfile.write("\n\n")

while True:
    line2 = infile2.readline()
    if len(line2) == 0: 
	break
    outfile.write(line2)

infile1.close()
infile2.close()
outfile.close()

print "Done!"
