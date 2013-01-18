#! /usr/bin/python

import sys
from gaspy import *
from radpy import *
from cfpylib.gasdyn.cea2_gas import *
from cfpylib.util.YvX import *
from time import time
from math import pi

def main():
    gas_input_file = sys.argv[1]
    rad_input_file = sys.argv[2]
    
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
    # Argon gas
    reactants["Ar"] = 1.0
    print reactants

    # calculate equilibrium gas state
    p = 1.0e5; T = 1.0e4
    Q = Gas_data(gm)
    Q.p = p
    for itm in range(ntm): Q.T[itm] = T
    for isp,sp in enumerate(species): Q.massf[isp] = get_species_composition(sp,reactants)
    gm.eval_thermo_state_pT(Q)
    cea = Gas( reactants, onlyList=reactants.keys(), with_ions=True, trace=1.0e-30 )
    cea.set_pT(p,T,transProps=False)
    #over-write provided initial mass-fractions
    #print "f_eq = ", eps.eq_massf
    for isp,sp in enumerate(species):
        Q.massf[isp] = get_species_composition(sp,cea.species)
    scale_mass_fractions( Q.massf )
    gm.eval_thermo_state_pT(Q)
    Q.print_values(False)

    print "temperature = ", Q.T[0]
    for isp in range(nsp):
	print "%s_partial_density = %e" % ( species[isp], Q.massf[isp]*Q.rho )
    print "\n"

    # now the tangent-slab calculation
    print "Performing a tangent slab transport calculation..."
    nslabs = 100
    s0 = 0.0; T0 = 0.0
    s1 = 0.2; T1 = 0.0
    ds = (s1-s0)/nslabs
    TS = TS_data(rsm,nslabs,T0,T1)
    divq_list = []; s_list = []
    for i in range(nslabs):
    	divq_list.append(new_doublep())
	s = s0+0.5*ds+ds*i
	s_list.append(s)
	print "creating new radpoint at s=%e" % s
    	TS.set_rad_point(i,Q,divq_list[-1],s,ds)
	j_total = TS.get_rpoint_pointer(i).X_.integrate_emission_spectra()
	print "divq(OT) = %e W/m3" % ( j_total * 4 * pi )
    q_total = TS.solve_for_divq_OT()
    print "Total radiative flux from a %f cm slab: q_total (OT) = %0.4f W/cm**2" % ( (s1-s0)*100, q_total*1.0e-4 )
    q_total = TS.quick_solve_for_divq()
    print "Total radiative flux from a %f cm slab: q_total (quick) = %0.4f W/cm**2" % ( (s1-s0)*100, q_total*1.0e-4 )
    q_total = TS.exact_solve_for_divq()
    print "Total radiative flux from a %f cm slab: q_total (exact) = %0.4f W/cm**2" % ( (s1-s0)*100, q_total*1.0e-4 )


    # write divq's to file
    ofile = open("divergence-profile.txt","w")
    ofile.write("# Column 1: Location, s (m)\n")
    ofile.write("# Column 2: Radiative divergence, divq (W/m3)\n")
    for i in range(nslabs):
	ofile.write("%e \t %e\n" % ( s_list[i], doublep_value(divq_list[i]) ) )
	delete_doublep(divq_list[i])
    ofile.close()    

    # compute average absorption coefficient
    B_vec = vectorSB()
    create_spectral_bin_vector( TS.get_rpoint_pointer(0).X_.kappa_nu, OPACITY_BINNING, 1, B_vec )
    X = BinnedCoeffSpectra( TS.get_rpoint_pointer(0).X_, B_vec )
    print "kappa_av = ", X.kappa_bin[0]

    # try to use YvX to plot the data
    TS.F_.write_to_file("flux_spectra.txt")
    qvL = YvX("flux_spectra.txt",0,1,False)
    # reverse the data and calculate the logarithm of j
    x_list = list(qvL.x_array); x_list.reverse(); qvL.x_array = array(x_list)
    y_list = list(qvL.y_array); y_list.reverse(); qvL.y_array = array(y_list)
    for i,y in enumerate(qvL.y_array):
        qvL.y_array[i] = log(y)
    qvL.plot_data(xlabel="wavelength, lambda (nm)", ylabel="Spectral flux, q (W/m**2-sr-m)", show_plot=True, include_integral=False)
    del TS

    del rsm, gm
    
    print "done."
    
if __name__ == '__main__':
    main()
