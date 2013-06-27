#! /usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Wed Apr 17 15:53:11 2013
Last update on Thurs 27 June, 2013 
    + Added electric dipole moment for CH
	+ Additional comments

authors: Nikhil Banerji <nikhil.banerji@epfl.ch>
         Jeremy Mora-Monteros <jeremy.mora-monteros@epfl.ch>

'''

import sys, math, os
from getopt import getopt, GetoptError
from calc_vars import *

shortOptions=""
longOptions=["help", "species=", "T=", "usr=", "reactype="]

def printUsage():
    print ""
    print "Usage:"
    print "CIC.py [--help] [--species = collision integral pair]"
    print "       [--T = temperature range vector]"
    print "       [--usr = account username to add to pathway to search repository for species files]"
    print "       [--reactype = reaction type - atom-atom (aa) or molecule-molecule(mm) or atom-molecule(am)]"
    print "e.g. ./cluster_CIC.py --species='N2 C2' --T='2000 5000 8000 11000 15000 20000 26000 32000' --usr='nbanerji' --reactype='mm'"
    return

# Lennard Jones/Stockmayer model used - only works with neutral particle interactions
# Not an issue as electron-neutral particle interactions do not contribute
# significantly to overall transport properties
# Also electrons & polyatomic species don't exist together in high concentrations

##---------------------------------------------------------------------##
# Read depth of intermolecular potential minimum (epsilon(Joules)) from 
# relevant species .lua files found in Eilmer3 repository
def read_eps(species):
    f1 = species[0] + ".lua"
    f2 = species[1] + ".lua"
    file1 = open(f1,"r")
    lines1 = file1.readlines()
    for i in range(0, len(lines1)):
        line1 = lines1[i]
        if species[0] + ".eps0" in line1:
            val1 = lines1[i+1].split("=")[-1]
            eps1 = float(val1.split(",")[0])
    file2 = open(f2,"r")
    lines2 = file2.readlines()
    for j in range(0,len(lines2)):
        line2 = lines2[j]
        if species[1] + ".eps0" in line2:
            val2 = lines2[j+1].split("=")[-1]
            eps2 = float(val2.split(",")[0])
    eps = [eps1, eps2]
    return eps  

def calc_eps_k(eps): # Calculate epsilon/kb
    k_b = [1.3806503e-23] # Boltzmann coefficient
    ek_ab = [(eps[0]*eps[1])**0.5]
    ek_ab = ek_ab[0]    
    e_k = ek_ab/k_b[0]
    #print "eps/kb = %0.3f" %(e_k)
    return e_k, ek_ab
        
def calc_T_r(T, e_k): # Calculate T*
    T_r = [x/e_k for x in T]
    #print "T* = " 
    #print " %0.2f "*len(T) % tuple(T_r) # reduce precision
    return T_r
    
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
    print "Calculating collision integrals..."
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
    print "species list:", species
    T = uoDict.get("--T","").split()
    for i in range(0,len(T)):
        T[i] = float(T[i])
    print "Temperature range:", T
    usr = uoDict.get("--usr", "") 
    reaction = uoDict.get("--reactype", "")
    if reaction == "aa":
        print "reaction type = atom-atom"
    elif reaction == "am":
        print "reaction type = atom-molecule"
    elif reaction == "mm":
        print "reaction type = molecule-molecule"
    
    # Now change the directory from the current directory
    os.chdir("/home/"+usr+"/cfcfd3/lib/gas/species/")
    
    # Read eps and sigma values from .lua species files
    # Calculate basic values
    eps = read_eps(species)
    e_k = calc_eps_k(eps)[0]
    ek_ab = calc_eps_k(eps)[1]
    T_r = calc_T_r(T, e_k)
    sigma = read_sigma(species)
    sigma_ab = calc_sigma_ab(sigma)
    print ek_ab
    
    
##--------------------------- Integration timestep values-----------#
# Optima calculated by Romain Savajano, EPFL, for reference please ask for 
# MATLAB files from IAG, EPFL.

    r_step = 1e-3   # optimum = 1e-3
    r_max = 1e3     # optimum = 1e3
    b_min = 0       # optimum = 0
    b_step = 0.1    # optimum = 0.1
    b_max = 2.5     # optimum = 2.5
    g_min = 0.1     # optimum = 0.1
    g_step = 2      # optimum = 2
    g_max = 50.1    # optimum = 50.1 guarantees an integer number of steps
    
##------------------------------ Calculate collision integrals ------#

# Calculation of delta for Stockmayer Potential
# To add further polar species, include in sp_polar and mu
    mu_m = [0]*len(species)
    # As no value was found for mu for C2H (polar), it can be set to zero, as it is "relatively small"
    # according to Tucker (1974)
    sp_polar = ['HCN','HCO','CN','CO','NO','CH'] 
    mu = {'HCN':2.985,'HCO':1.35,'CN':1.45,'CO':0.112,'NO':0.0672,'CH':0.88} # Values from Krieger et al.(1951), Scarl et al.(1974) in Debye
    if species[0] in sp_polar:
        polar = species[0]
        if mu.has_key(polar):
            mu_m[0] = mu[polar]
    if species[1] in sp_polar:
        polar = species[1]
        if mu.has_key(polar):
            mu_m[1] = mu[polar]
    # 1e-10 to convert to esu-angstrom; 1e10 for m to angstrom.
    delta = (1e-10)**2*(mu_m[0]*mu_m[1])/(2*ek_ab*(sigma_ab[0]*1e10)**3)
    print "delta =", delta
    
    # progress bar - displayed in outfile
    toolbar_width = 2*len(T_r)
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line '['
    
    # Following calculation done using functions described in 
    # calc_vars -> l,s = [1,1] for diffusion CI
    omega11 = [0]*(len(T_r))
    omega11_r = [0]*(len(T_r))    
    for i in range(0,len(T_r)):
        [l,s] = [1,1]
        omega11[i] = Omegals_r(l,s,T_r[i],g_min,g_step,g_max,b_min,b_step,b_max,r_step,r_max,delta)        
    # At higher temperatures, collision integral is badly approximated using
    # the method used. Therefore Park et al. and Kim et al. proposed
    # the following corrections.            
        if reaction == "aa" or reaction == "am": # Use Park et al. correlation
            omega11_r[i] = 0.9*omega11[i]*(T[i]/2000)**-0.25
        elif reaction == "mm": # Use Kim et al. correlation
            omega11_r[i] = omega11[i]*(0.5+0.45*(T[i]/2000)**-0.25)
        else:
            omega11_r = omega11
        sys.stdout.write("-")
        sys.stdout.flush()
    # l,s = [2,2] for viscosity CI
    omega22 = [0]*(len(T_r))
    omega22_r = [0]*(len(T_r))    
    for j in range(0,len(T_r)):
        [l,s] = [2,2]
        omega22[j] = Omegals_r(l,s,T_r[j],g_min,g_step,g_max,b_min,b_step,b_max,r_step,r_max,delta)
            
        if reaction == "aa" or reaction == "am":
            omega22_r[j] = 0.9*omega22[j]*(T[j]/2000)**-0.25
        elif reaction == "mm":
            omega22_r[j] = omega22[j]*(0.5+0.45*(T[j]/2000)**-0.25)
        else:
            omega22_r = omega22
            
        sys.stdout.write("-")
        sys.stdout.flush()
    sys.stdout.write("\n")
    
    print "omega(1,1)* w/o correction ="
    print " %0.4f "*len(omega11) % tuple(omega11) # reduce precision
    print "omega(2,2)* w/o correction ="
    print " %0.4f "*len(omega22) % tuple(omega22) # reduce precision
    print
    print "omega(1,1)* ="
    print " %0.4f "*len(omega11_r) % tuple(omega11_r) # reduce precision
    print "omega(2,2)* ="
    print " %0.4f "*len(omega22_r) % tuple(omega22_r) # reduce precision

# To launch simulations for each temperature separately 
# (simulate mpi - quicker simulations on cluster, see shell script), 
# post-proc done using CIC_post.py
# If not and simulation launched as in example at the top of the page,
# please un-comment remainder of the code.

'''	
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
    # Can only write in scratch dir on cluster    
    os.chdir("/scratch/"+usr+"/collision_integrals")  
    #os.chdir("/home/"+usr+"/Documents/collision_integrals/CIC")
    filename = uoDict.get("--species", "")
    # Name CI file in eilmer3 compatible format
    if species[0]+" " in filename:
        filename = species[0]+"-"+species[1]+".lua"
    # Write CI file
    output = open(filename,'w')
    output.write("CI = {\n")
    output.write("  i = ''%s'',\n" % species[0])
    output.write("  j = ''%s'',\n" % species[1])
    output.write("  reference = 'EPFL, 2013 - Collision Integral Calculator - Unpublished',\n")
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
'''    
##------------------------------------------------------------------------##

if __name__ == '__main__':
    main()
