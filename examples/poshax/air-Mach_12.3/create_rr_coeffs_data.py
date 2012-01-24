#!/usr/bin/env python
# create_rr_coeffs_data.py

import sys
from libgas2 import *

def printUsage():
    print "create_rr_coeffs_data.py"
    print "Creates data for the reaction rate coefficients over a specified"
    print "temperature range for a given reaction scheme."
    print "Usage:"
    print "    create_rr_coeffs_data.py <themo_file> <reactions_file> T_low T_high <outfile>"
    print "  where"
    print "    <thermo_file>    : a thermodynamic inupt file"
    print "    <reactions_file> : a reaction scheme input file"
    print "    T_low            : lower limit of temperature range"
    print "    T_high           : upper limit of temperature range"
    print "    <outfile>        : file for output data file"
    sys.exit(1)

def main():
    if len(sys.argv) != 6:
        printUsage()
        
    thermo_file = sys.argv[1]
    reac_file = sys.argv[2]
    T_low = float(sys.argv[3])
    T_high = float(sys.argv[4])
    fp = open(sys.argv[5], 'w')

    if T_low > T_high:
        print "T_low= ", T_low, " is greater than T_high= ", T_high
        printUsage()

    dT = (T_high - T_low) / 1000.0

    # Initialise the gas model and reaction scheme
    set_type_of_gas("perf_gas_mix", thermo_file)
    set_reaction_scheme( "test scheme",  reac_file )

    Q = gas_data()
    r = get_reaction_scheme_pointer()
    reac_list = []
    for i in range(r.get_number_of_reactions()):
        reac_list.append( r.get_reaction_pointer(i) )
       
    T = T_low
    while T <= T_high:
        Q.T = T
        print "T= ", T
        fp.write("%20.12e" % T)
        for ir in range(r.get_number_of_reactions()):
            reac = reac_list[ir]
            reac.fill_in_forward_coeff( Q )
            reac.fill_in_backward_coeff( Q )
            k_f = reac.k_f()
            k_b = reac.k_b()
            fp.write(" %20.12e %20.12e" % (k_f, k_b) )
        fp.write("\n")
        T += dT
    

if __name__ == '__main__':
    main()

        

    
