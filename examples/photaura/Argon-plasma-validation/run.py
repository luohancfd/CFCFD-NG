#! /usr/bin/python

import sys
from gaspy import *
from radpy import *
from cea2_gas import *

# make the gas-model
species = ["Ar","Ar_plus","e_minus"]
gm = create_gas_file( "one temperature gas", species, "gas-model.lua" )
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()

# create the radiation model
rsm = create_radiation_spectral_model( "rad-model.lua" )

massfs = [ 1.0, 0.0, 0.0 ]
    
# gas pressure
p = PC_P_atm

# result lists
T_list = []
j_list = []

for T in range(5000,15001,1000):
    # setup the reactants list
    reactants = make_reactants_dictionary( species )
    for sp in species:
        _sp = sp.replace("_plus","+").replace("_minus","-")
        reactants[_sp] = massfs[gm.get_isp_from_species_name(sp)]
    
    massf_sum = 0.0
    for massf in reactants.values():
        massf_sum += massf
    for sp in reactants.keys():
        reactants[sp] /= massf_sum 
    print reactants

    # solve the pT equilibrium radiation problem
    Q = Gas_data(gm)
    Q.p = p
    for itm in range(ntm): Q.T[itm] = T
    for isp,sp in enumerate(species): Q.massf[isp] = get_species_composition(sp,reactants)
    cea = Gas( reactants, onlyList=reactants.keys(), with_ions=get_with_ions_flag(species), trace=1.0e-20 )
    cea.set_pT(p,T,transProps=False)
    # print "filtered species composition: ", cea.species
    #over-write provided initial mass-fractions
    for isp,sp in enumerate(species):
        Q.massf[isp] = get_species_composition(sp,cea.species)
            
    gm.eval_thermo_state_pT(Q)
    # print "computed equlibrium state with filtered species: "
    # Q.print_values(False)

    if "e_minus" in species:
        N_elecs = ( Q.massf[gm.get_isp_from_species_name("e_minus")]*Q.rho/RC_m_SI*1.0e-6 )
        N_total = ( ( Q.p - Q.p_e ) / RC_R_u / Q.T[0] + Q.p_e / RC_R_u / Q.T[-1] ) * RC_Na * 1.0e-6
        print "\ngas temperature = %e K" % T
        print "electron number density = %e cm-3" % ( N_elecs )
        print "total number density = %e cm-3" % N_total
        print "ionization fraction = %e" % ( N_elecs / N_total )

    for sp in species:
        isp = gm.get_isp_from_species_name(sp)
        N = ( Q.massf[isp]*Q.rho/(gm.molecular_weight(isp)/RC_Na)*1.0e-6 )
        print "%s number density = %e cm-3" % ( sp, N )

    # perform spectra calculation
    j_total = rsm.radiative_integrated_emission_for_gas_state( Q, True )
    j_list.append( j_total )
    T_list.append (T )

ofile = open("j_versus_T.txt","w")
for i,T in enumerate(T_list):
    j = j_list[i]
    ofile.write("%e %e\n" % ( T, j ) )
ofile.close()

print "Done."
