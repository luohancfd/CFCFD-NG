#! /usr/bin/python

import sys
from gaspy import *
from librad2 import *
try:
    from cfpylib.gasdyn.cea2_gas import *
    cea2_gas = True
except:
    print "WARNING: The cea2_gas module is not functional."
    print "         Equilibrium calculation disabled."
    cea2_gas = False
    
from time import time
try:
    from cfpylib.util.YvX import *
except ImportError:
    print "Problem loading YvX module."
    print "Make sure cfpylib is in your installation directory."
    sys.exit()

def create_gas_file(model, species, fname):
    # borrowed from "e3_gas.py"
    fp = open("gas-tmp.inp", "w")
    fp.write("model = '%s'\n" % model)
    fp.write("species = {")
    for spec in species:
        fp.write("'%s', " % spec)
    fp.write("}\n")
    fp.close()

    cmd = "%s %s %s" % ("gasfile", "gas-tmp.inp", fname)
    os.system(cmd)
    cmd = "rm gas-tmp.inp"
    os.system(cmd)

    return

EMISSION_CALC = True
TS_CALC = True

def main():
    if len(sys.argv) > 1:
        gas_input_file = sys.argv[1]
	rad_input_file = sys.argv[2]
    else:
        gas_input_file = "11sp-air-gm.lua"
	rad_input_file = "air-radiators.lua"
	
    species_list = [ "N2", "N2_plus", "O2", "O2_plus", "NO", "NO_plus", "O", "O_plus", "N", "N_plus", "e_minus" ]
    create_gas_file( "two temperature gas", species_list, gas_input_file )
	
    print "Setting up a two temperature gas mix from thermo file", gas_input_file 
    gm = create_gas_model(gas_input_file)

    print "Setting up a managed radiation model from input file", rad_input_file
    rsm = create_radiation_spectral_model( rad_input_file )
    print "\n"
    
    # species  =  N2   O2   NO   N    O 
    # index  =    0    1    2    3    4
    
    nsp = gm.get_number_of_species()
    ntm = gm.get_number_of_modes()
    Q = Gas_data(nsp,ntm)
    if cea2_gas:
    	# Compute equilibrium composition via CEA 
        species  =  []
        for isp in range(nsp):
    	    species.append(gm.species_name(isp))
    	    species[-1] = species[-1].replace("_plus","+")
    	    species[-1] = species[-1].replace("_minus","-")
        f_inf = [ 0.0 ] * nsp
        f_inf[gm.get_isp_from_species_name('N2')] = 0.767; f_inf[gm.get_isp_from_species_name('O2')] = 0.233
        p_inf = 40.0; T_inf = 296.0; u_inf = 10.254e3
        Q.p = p_inf
        for itm in range(ntm): Q.T[itm] = T_inf
        for isp in range(nsp): Q.massf[isp] = f_inf[isp]
        gm.eval_thermo_state_pT(Q)
        cea = Gas('mix', species, f_inf, use_out_file=True)
        cea.set_from_pAndT(p_inf,T_inf)
        eps = cea.shock_process( u_inf )
        #over-write provided initial mass-fractions
        #print "f_eq = ", eps.eq_massf
        Q.rho = eps.rho
        for itm in range(ntm): Q.T[itm] = eps.T
        for isp in range(nsp): Q.massf[isp] = eps.eq_massf[isp]
	del cea, eps
    else:
	# Use precomputed equilibrium composition as cea2_gas is not functioning
        Q.rho = 7.2917e-03
        for itm in range(ntm): Q.T[itm] = 10441.0
        eq_massf = [ 1.43390e-03, 6.50260e-05, 1.2349e-06, 4.6649e-07, 6.1227e-05,
                     1.10000e-04, 2.24290e-01, 8.6206e-03, 7.2153e-01, 4.3886e-02,
                     2.01780e-06 ]
 	for isp in range(nsp): Q.massf[isp] = eq_massf[isp]
    gm.eval_thermo_state_rhoT(Q)
    Q.print_values(False)

    # print "electron number density = %e cm-3" % ( Q.massf[gm.get_isp_from_species_name("e-")]*Q.rho/RC_m_SI*1.0e-6 )
    if EMISSION_CALC:
	# perform radiation calculation
	n_evals = 1
	print "\nPerforming %d emission coefficient calculations..." % n_evals
	t0 = time()
	for i in range(n_evals): j_total = rsm.radiative_integrated_emission_for_gas_state( Q, True )
	t1 = time()
	print "Wall time = %f seconds" % ( t1-t0 ) 
	print "j_total = %0.2e W/cm**3\n" % ( j_total * 1.0e-6 )
	
	print "Calculating emission and absorption spectra for gas state"
	X = CoeffSpectra()
	rsm.radiative_spectra_for_gas_state( Q, X )
	j_total = X.write_to_file( "coefficient_spectra.txt" );
	print "Spectrally resolved integrated emission: j_total = %0.2f W/cm**3" % ( j_total * 1.0e-6 )
	j_total = rsm.radiative_integrated_emission_for_gas_state( Q, False )
	print "Spectrally unresolved integrated emission: j_total = %0.2f W/cm**3" % ( j_total * 1.0e-6 )
	
	# try to use YvX to plot the data
	JvL = YvX("coefficient_spectra.txt",0,1,False)
	# reverse the data and calculate the logarithm of j
	x_list = list(JvL.x_array); x_list.reverse(); JvL.x_array = array(x_list)
	y_list = list(JvL.y_array); y_list.reverse(); JvL.y_array = array(y_list)
	for i,y in enumerate(JvL.y_array):
	    JvL.y_array[i] = log(y)
	JvL.plot_data(xlabel="wavelength, lambda (nm)", ylabel="Spectral emission coefficient, j (W/m**3-sr-m)", show_plot=True, include_integral=False)
    
	# show how certain properties can be outputed to file
	rsm.write_line_widths_to_file(Q)

    if TS_CALC:
	# now demonstrate use of the LOS/TS classes to perform a tangent-slab calculation
	print "Performing a tangent slab transport calculation..."
	nslabs = 1
	s0 = 0.0; T0 = 0.0
	s1 = 5.0e-2; T1 = 1000.0
	ds = (s1-s0)/nslabs
	TS = TS_data(rsm,nslabs,T0,T1)
	for islab in range(nslabs):
	    divq = new_doublep()
	    TS.set_rad_point(islab,Q,divq,s0+0.5*ds+islab*ds,ds)
	q_total = TS.quick_solve_for_divq()
	print "Total radiative flux from a %f cm slab: q_total = %0.2f W/cm**2" % ( (s1-s0)*100, q_total*1.0e-4 )
	print "divq = %0.2f W/m**3" % doublep_value(divq)
	delete_doublep(divq)
	
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

    del gm, rsm, Q
    
    print "done."
    
if __name__ == '__main__':
    main()
