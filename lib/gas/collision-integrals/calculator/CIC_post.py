#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:54:24 2013

@author: Jeremy Mora-Monteros
"""

import sys, math, os
from getopt import getopt, GetoptError
import numpy as np

shortOptions=""
longOptions=["help", "species=", "T=", "usr="]

def printUsage():
    print ""
    print "Usage:"
    print "CIC_post.py [--help] [--species = collision integral pair]"
    print "       [--T = temperature range vector]"
    print "       [--usr = account username to add to pathway to search repository for species files]"
    print "e.g. ./CIC_post.py --species='N2 C2' --T='2000 4000 8000 16000 32000' --usr='nbanerji'"
    return
    
##--------------------------reading outfiles function-------------------#

def read_omega(species,T):
    f = [0]*len(T) 
    omega11 = [0]*len(T)
    omega22 = [0]*len(T)
    for i in range(0, len(T)):
        f[i] = "outfile_" + species[0] + "-" + species[1] + "_" + T[i]
    for i in range(0, len(f)):
        file1 = open(f[i],"r")
        lines1 = file1.readlines()
        line1 = lines1[i]
        for j in range(0, len(lines1)):
            line1 = lines1[j]
            if "omega(1,1)* =" in line1:
                val1 = lines1[j+1].split("=")[-1]
                omega1 = float(val1.split(",")[0])
            elif "omega(2,2)* =" in line1:
                val2 = lines1[j+1].split("=")[-1]
                omega2 = float(val2.split(",")[0])
        omega11[i] = omega1
        omega22[i] = omega2
        file1.close()
    return omega11, omega22
    
## Calculate sigma for both species

def read_sigma(species):
    f1 = species[0] + ".lua"
    f2 = species[1] + ".lua"
    file1 = open(f1,"r")
    lines1 = file1.readlines()
    for i in range(0, len(lines1)):
        line1 = lines1[i]
        if species[0] + ".sigma" in line1:
            val1 = lines1[i+1].split("=")[-1]
            sigma1 = float(val1.split(",")[0])
    file2 = open(f2,"r")
    lines2 = file2.readlines()
    for j in range(0,len(lines2)):
        line2 = lines2[j]
        if species[1] + ".sigma" in line2:
            val2 = lines2[j+1].split("=")[-1]
            sigma2 = float(val2.split(",")[0])
    sigma = [sigma1, sigma2]
    return sigma
            
## Sigma AB
def calc_sigma_ab(sigma):
    sigma_ab=[(sigma[0]+sigma[1])/2]
    return sigma_ab
    
##------------------------------- main function ------------------------#

def main():  
    # Read inputs from command line
    try:
        userOptions = getopt(sys.argv[1:],shortOptions, longOptions)
    except GetoptError, e:
        print "One (or more) of your command-line options was no good."
        print "     ", e
        printUsage()
        sys.exit(1)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(0)
        #
    species = uoDict.get("--species", "").split()
    print "Create .lua file for the collision integral pair", species
    
    T = uoDict.get("--T","").split()
    print "Temperature range:", T
    
    usr = uoDict.get("--usr", "") 
    
    od = "outfiles_" + species[0] + "-" + species[1]
    print "Outfile directory:", od
    
    # Now change the directory from the current directory
    os.chdir("./"+od)
    
    # Read Omega* from outfiles
    omega11,omega22 = read_omega(species,T)
    
# Convert T values to floats
    for i in range(0,len(T)):
        T[i] = float(T[i])

    os.chdir("/home/"+usr+"/cfcfd3/lib/gas/species/")

    # Read sigma values from .lua species files
    sigma = read_sigma(species)
    sigma_ab = calc_sigma_ab(sigma)
    print sigma_ab

# Calculate Gupta-Yos curve fit coefficients for Eilmer3 input
    Y11 = [0]*(len(T))
    Y22 = [0]*(len(T))
    T_log = [0]*(len(T))     
    for k in range(0,len(T)):
        Y11[k] = np.log(omega11[k]*math.pi*sigma_ab[0]**2)
        Y22[k] = np.log(omega22[k]*math.pi*sigma_ab[0]**2)
        T_log[k] = np.log(T[k])
    P11 = np.polyfit(T_log,Y11,3)
    P22 = np.polyfit(T_log,Y22,3)
    
    print ""        
    print "omega11_pi =", P11
    print ""
    print "omega22_pi =", P22
##------------------------Write output file for Eilmer3 -------------#

    # Now change the directory to where you want the CI files to be written
    # Can only write in scratch dir on cluster at EPFL    
    os.chdir("/scratch/"+usr+"/collision_integrals")  
    #os.chdir("../")
    filename = uoDict.get("--species", "")
    # Name CI file in eilmer3 compatible format
    if species[0]+" " in filename:
        filename = species[0]+"-"+species[1]+".lua"
    # Write CI file
    output = open(filename,'w')
    output.write("CI = {\n")
    output.write("  i = ''%s'',\n" % species[0])
    output.write("  j = ''%s'',\n" % species[1])
    output.write("  reference = 'EPFL, 2013 - Collision integral calculator - Unpublished',\n")
    output.write("  model = 'GuptaYos curve fits',\n")
    output.write("  parameters = {\n")
    output.write("    {\n")
    output.write("      T_low = %s,\n" % T[0])   
    output.write("      T_high = %s,\n" % T[len(T)-1])
    output.write("      Pi_Omega_11 = {      %.5f,       %.5f,       %.5f,       %.5f, },\n" %(P11[0], P11[1], P11[2], P11[3]))
    output.write("      Pi_Omega_22 = {      %.5f,       %.5f,       %.5f,       %.5f, },\n" %(P22[0], P22[1], P22[2], P22[3]))
    output.write("      D           = {      %.5f,       %.5f,       %.5f,       %.5f, },\n" %(0, 0, 0, 0))
    output.write("    },\n")
    output.write("  },\n")
    output.write("},\n")
    
##------------------------------------------------------------------------##

if __name__ == '__main__':
    main()
