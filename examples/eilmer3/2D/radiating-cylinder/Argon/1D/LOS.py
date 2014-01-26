#! /usr/bin/python

import sys
from gaspy import *
from radpy import *
from cfpylib.gasdyn.cea2_gas import *
from cfpylib.util.YvX import *
from time import time
from math import pi

def main():
    if len(sys.argv) > 1:
        gas_input_file = sys.argv[1]
        rad_input_file = sys.argv[2]
    else:
        gas_input_file = "gas-model.lua"
        rad_input_file = "rad-model.lua"
        create_gas_file( "one temperature gas", ["Ar", "Ar_plus", "e_minus" ], gas_input_file )
	
    print "Setting up a thermally perfect gas mix from thermo file", gas_input_file 
    gm = create_gas_model(gas_input_file)

    print "Setting up a managed radiation model from input file", rad_input_file 
    rsm = create_radiation_spectral_model( rad_input_file )
    
    # may want to compute equilibrium composition via CEA
    species  =  []
    nsp = gm.get_number_of_species()
    ntm = gm.get_number_of_modes()
    for isp in range(nsp):
        species.append(gm.species_name(isp))
    reactants = make_reactants_dictionary( species )
    # assume Argon
    reactants["Ar"] = 1.0
    massf_sum = 0.0
    for massf in reactants.values():
        massf_sum += massf
    for sp in reactants.keys():
        reactants[sp] /= massf_sum 
    print reactants

    cea = Gas( reactants, onlyList=reactants.keys(), with_ions=True, trace=1.0e-10 )

    # gas slab 1
    L = 10.0e-2
    p_inf = 1.0e5; T_inf = 10.0e3
    cea.set_pT(p_inf,T_inf)
    #print "f_eq = ", eps.eq_massf
    Q = Gas_data(gm)
    Q.rho = cea.rho
    for itm in range(ntm): Q.T[itm] = cea.T
    for isp,sp in enumerate(species):
        Q.massf[isp] = get_species_composition(sp,cea.species)
    gm.eval_thermo_state_rhoT(Q)
    Q.print_values(False)
     
    # do the LOS calculations
    Q_rE_rad = new_doublep()
    TS = TS_data( rsm, 1 )
    TS.set_rad_point(0,Q,Q_rE_rad,L*0.5,L)
    j_total = TS.get_rpoint_pointer(0).X_.write_to_file("coefficient_spectra.txt")
    print "j_total inside slab = ", j_total
    q_total = TS.exact_solve_for_divq()
    q_total = TS.F_.write_to_file("flux_spectra.txt")
    print "q_total at end of slab = ", q_total

    

    print "done."
    
if __name__ == '__main__':
    main()

