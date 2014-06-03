#! /usr/bin/python

import sys
import os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
sys.path.append("") # so that we can find user's scripts in current directory
from gaspy import *
from time import time

if not os.path.exists('gasmodels'):
    os.makedirs('gasmodels')
if not os.path.exists('libgas_results'):
    os.makedirs('libgas_results')


def main():
    if len(sys.argv)!=2:
        print "usage: noneq-species-test.py <species-name>" 
    sp = sys.argv[1]
    print "Setting up a four temperature gas model for species: %s" % ( sp )
    gm_path = "gasmodels/" + sp+"-noneq-gm.lua"
    create_gas_file( "four temperature gas", [sp], gm_path )
    gm = create_gas_model(gm_path)
    nsp = gm.get_number_of_species()
    ntm = gm.get_number_of_modes()
    # prepare to loop over the following temperature range
    p = 1.0e5; T_min = 50.0; T_max = 50000.0; dT = 50.0
    Q = Gas_data(nsp,ntm)
    Q.p = p; T = T_min
    # prepare output files
    columns  = "# Column 1: temperature, T (K)\n"
    columns += "# Column 2: specific heat at constant pressure, Cp (J/mol-k)\n"
    columns += "# Column 3: specific heat at constant volume, Cv (J/mol-k)\n"
    columns += "# Column 4: ratio of specific heats, gamma (ND)\n"
    ofile_name = sp + ".txt"
    ofile = open( "libgas_results/" + ofile_name, "w" )
    ofile.write("# Filename: %s\n" % ofile_name)
    ofile.write(columns)
    t0 = time()
    while T<=T_max:
        # print "computing properties for: T = ", T
        for itm in range(ntm): Q.T[itm] = T
        Q.massf[0] = 1.0
        gm.eval_thermo_state_pT(Q)
        Cp = gm.Cp(Q) * gm.molecular_weight(0)
        Cv = gm.Cv(Q) * gm.molecular_weight(0)
        gamma = gm.gamma(Q)
        ofile.write("%0.4e \t %0.4e \t %0.4e \t %0.4e\n"\
                     % ( Q.T[-1], Cp, Cv, gamma ) )
        T += dT
    t1 = time()
    ofile.close()
    call_gas_model_deconstructor(gm)     
    print "\nElapsed time for species %s: %0.2e seconds" % ( sp, t1-t0 )
    print "\ndone."
    
if __name__ == '__main__':
    main()
