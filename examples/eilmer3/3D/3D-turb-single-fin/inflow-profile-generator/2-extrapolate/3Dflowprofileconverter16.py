#!/usr/bin/env python
# 3Dinfileconverter.py
# Accepts 2D BL Profile, extrapolates into 3D for use as an input file for StaticProfBC()
# Based on compute_viscous_data_simple.py by Peter Jacobs and Wilson Chan
# Created 17 August 2014 by Samuel Stennett

#EXECUTE USING python 3Dinfileconverter.py OTHERWISE WILL NOT WORK

ncelly = 25     # Number of cells in y direction in 3D extrapolation
		# Compensate for block splitting! 

ncellz = 25	# Number of cells in z direction, before break outfile (to compensate for block splitting)

print "Begin 3D Converter"

infile = open("./interpolated_profile_ghost1.dat", "r")
infile2 = open("./interpolated_profile_ghost2.dat","r")
line = infile.readline()  # read away the first (comment) line containing notes on format
line2 = infile2.readline()
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile1 = open("3Dinfilelower1.dat", "w")
outfile2 = open("3Dinfilelower2.dat", "w")
outfile3 = open("3Dinfileupper1.dat", "w")
outfile4 = open("3Dinfileupper2.dat", "w")
outfile1.write("# pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0] \n")
outfile2.write("# pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0] \n")
outfile3.write("# pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0] \n")
outfile4.write("# pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0] \n")

# Process every data line.
i=0
while i < ncellz:
    line = infile.readline().strip()
    line2 = infile2.readline().strip()
    if len(line) == 0: break
    if len(line2) == 0: break
    init2 = 0
    while init2 < ncelly:
        outfile1.write(line + "\n")
        outfile1.write(line2 + "\n")
        init2 += 1
    i += 1
i=0
while i < ncellz:
    line = infile.readline().strip()
    line2 = infile2.readline().strip()
    if len(line) == 0: break
    if len(line2) == 0: break
    init2 = 0
    while init2 < ncelly:
        outfile2.write(line + "\n")
        outfile2.write(line2 + "\n")
        init2 += 1
    i += 1
i=0
while i < ncellz:
    line = infile.readline().strip()
    line2 = infile2.readline().strip()
    if len(line) == 0: break
    if len(line2) == 0: break
    init2 = 0
    while init2 < ncelly:
        outfile3.write(line + "\n")
        outfile3.write(line2 + "\n")
        init2 += 1
    i += 1
i=0
while i < ncellz:
    line = infile.readline().strip()
    line2 = infile2.readline().strip()
    if len(line) == 0: break
    if len(line2) == 0: break
    init2 = 0
    while init2 < ncelly:
        outfile4.write(line + "\n")
        outfile4.write(line2 + "\n")
        init2 += 1
    i += 1
	
infile.close()
infile2.close()
outfile1.close()
outfile2.close()
outfile3.close()
outfile4.close()

print "Done."
