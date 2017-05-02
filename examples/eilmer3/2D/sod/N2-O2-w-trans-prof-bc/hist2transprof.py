#!/usr/bin/env python
#
# Author: Rowan J. Gollan
# Date: 2017-04-28
#
# This script is used to convert a file with two adjacent history cells
# into a transient profile boundary condition data over a slice of cells.

ncellsInProfile = 2
inpFile = "hist/sod.hist.b0000"
outFile = "sod-driver-transient-profile-bc-q.dat"

ofp = open(outFile, "w")
ofp.write("# data converted from hist file.\n")
ofp.write("# Variables list: x y z volume rho u v w p a mu k0 mu_t k_t S tke omega mf0 mf1 dt_chem e0 T0\n")
ofp.write("# 2\n")

ifp = open(inpFile, "r")
ifp.readline() # Throw away header line

# Now begin looping through file.
# We'll work in pairs of lines and then expand that 
# out for the full profile in the out file.

while True:
    rghtCell = ifp.readline()
    if not rghtCell:
        break
    rTks = rghtCell.split()
    leftCell = ifp.readline()
    lTks = leftCell.split()
    t = rTks[0]
    ofp.write("time= %20.12e\n" % t) #as in block_io.cxx
    for i in range(ncellsInProfile):
        ofp.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %
                  (rTks[4], rTks[5], rTks[6], rTks[7], rTks[8], rTks[9], rTks[10], rTks[11],
                   rTks[12], rTks[13], rTks[14], rTks[15], rTks[16], rTks[17], rTks[18], rTks[19], rTks[20],
                   rTks[21], rTks[22], rTks[23], rTks[24], rTks[25]))
        ofp.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %
                  (lTks[4], lTks[5], lTks[6], lTks[7], lTks[8], lTks[9], lTks[10], lTks[11],
                   lTks[12], lTks[13], lTks[14], lTks[15], lTks[16], lTks[17], lTks[18], lTks[19], lTks[20],
                   lTks[21], lTks[22], lTks[23], lTks[24], lTks[25]))

ifp.close()
ofp.close()
                   
        
        





