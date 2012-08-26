#! /usr/bin/python

import sys
import os
from gaspy import *
from cfpylib.gasdyn.cea2_gas import *
from time import time
from random import random

def create_gas_file_list():
    gas_input_files = { "5sp-EQ" : "5sp-EQ.lua", "5sp-NEQ" : "5sp-NEQ.lua", "11sp-NEQ" : "11sp-NEQ.lua" }
    species = {}
    species["5sp"] = [ "N2", "O2", "NO", "N", "O" ]
    species["11sp"] = [ "N2", "N2_plus", "O2", "O2_plus", "NO", "NO_plus", "N", "N_plus", "O", "O_plus", "e_minus" ]
    gas_models = {}
    gas_models["EQ"] = "thermally perfect gas"
    gas_models["NEQ"] = "two temperature gas"
    for gm_name in gas_input_files.keys():
        i = gm_name.find("-")
        gm_type = gas_models[gm_name[i+1:]]
        species_list = species[gm_name[:i]]
        create_gas_file(gm_type,species_list,gas_input_files[gm_name])
        clear_chemical_species_library()
    return gas_input_files

def make_reactants_dictionary( species_list ):
    nsp = len(species_list)
    reactants = dict()
    for sp in species_list:
        # replace names containing '_plus' with '+' 
        if ( sp.find("_plus")>=0 ): sp = sp[0:sp.find("_plus")] + "+"
        # replace names containing '_minus' with '-' 
        if ( sp.find("_minus")>=0 ): sp = sp[0:sp.find("_minus")] + "-"
    reactants.setdefault(sp,0.0)
    return reactants

def get_species_composition( sp, species_data ):
    # replace names containing '_plus' with '+' 
    if ( sp.find("_plus")>=0 ): sp = sp[0:sp.find("_plus")] + "+"
    # replace names containing '_minus' with '-' 
    if ( sp.find("_minus")>=0 ): sp = sp[0:sp.find("_minus")] + "-"
    if sp in species_data.keys():
        return species_data[sp]
    else:
        return 0.0

def main():
    gas_input_files = create_gas_file_list()
    for gm_type in gas_input_files.keys():
        print "Setting up a %s gas model from thermo file: %s" % ( gm_type, gas_input_files[gm_type] )
        gm = create_gas_model(gas_input_files[gm_type] )
        species  =  []
        nsp = gm.get_number_of_species()
        ntm = gm.get_number_of_modes()
        for isp in range(nsp):
            species.append(gm.species_name(isp))
        # setup CEA gas
        reactants = make_reactants_dictionary( species )
        reactants["N2"] = 0.767; reactants["O2"] = 0.233
        if "e_minus" in species: with_ions=True
        else: with_ions=False
        cea = Gas( reactants, onlyList=reactants.keys(), with_ions=with_ions, trace=1.0e-10 )
        # prepare to loop over the following temperature range
        p = 1.0e5; T_min = 200.0; T_max = 13400.0; dT = 1000.0
        Q = Gas_data(gm)
        Q.p = p; T = T_min
        # prepare output files
        columns  = "# Column 0: iteration\n# Column 1: T\n# Column 2: rho\n# Column 3: enthalpy\n# Column 4: energy\n"
        columns += "# Column 5: entropy\n# Column 6: Cp\n# Column 7: gamma\n# Column 8: a\n"
        gfile_name = gm_type + "_libgas_results.txt"
        gfile = open( gfile_name, "w" )
        gfile.write("# Filename: %s\n" % gfile_name)
        gfile.write(columns)
        cfile_name = gm_type + "_cea2_results.txt"
        cfile = open( cfile_name, "w" )
        cfile.write("# Filename: %s\n" % cfile)
        cfile.write(columns)
        i = 0
        while T<=T_max:
            print "computing properties for: T = ", T
            cea.set_pT(p,T)
            #over-write provided initial mass-fractions
            for itm in range(ntm): Q.T[itm] = T
            for isp,sp in enumerate(species):
                Q.massf[isp] = get_species_composition(sp,cea.species)
            gm.eval_thermo_state_pT(Q)
            # perturb temperature and perform rhoe eval to test newton iterations
            for itm in range(ntm): Q.T[itm] = T*(0.5+random())
            gm.eval_thermo_state_rhoe(Q)
            # write gas and cea2 thermo properties to file
            h = gm.mixture_enthalpy(Q)
            u = Q.e[0]
            s = gm.mixture_entropy(Q)
            Cp = gm.Cp(Q)
            gam = gm.gamma(Q)
            gfile.write("%d %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e\n"\
                         % ( i, Q.T[-1],       Q.rho,   h,       u,  s,       Cp,      gam,     Q.a ) )
            cfile.write("%d %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e\n"\
                         % ( i, T,       cea.rho, cea.h,   cea.u,   cea.s,   cea.cp,  cea.gam, cea.son ) )
            # print "libgas.h / cea2.h = ", u/cea.u
            T += dT; i+=1

        gfile.close()
        cfile.close()
        
        nevals = 1000
        print "Timing %d rhoe evals at the last thermo state...\n" % nevals
        t0 = time();
        for i in range(nevals):
            # perturb temperature so we don't have too good a guess
            for itm in range(ntm): Q.T[itm] = (T-dT)*(0.5+random())
            gm.eval_thermo_state_rhoe(Q) 
        t1 = time(); print "Wall time = %f seconds" % ( t1-t0 ) 
            
        del gm
        clear_chemical_species_library()
   
    print "done."
    
if __name__ == '__main__':
    main()
