#! /usr/bin/python

import sys
from gaspy import *
from radpy import *
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
        reactants = make_reactants_dictionary( species_list )
        reactants['N2'] = 0.767; reactants['O2'] = 0.233
        p_inf = 40.0; T_inf = 296.0; u_inf = 10.254e3
        cea = Gas( reactants, onlyList=reactants.keys(), with_ions=True, trace=1.0e-10 )
        cea.set_pT(p_inf,T_inf)
        cea.shock_process( u_inf )
        Q.rho = cea.rho
        for itm in range(ntm):
            Q.T[itm] = cea.T
        for isp,sp in enumerate(species_list):
            Q.massf[isp] = get_species_composition(sp,cea.species)
        del cea
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
        j_total = X.write_to_file( "coefficient_spectra.txt", FREQUENCY );
        print "Spectrally resolved integrated emission: j_total = %0.2f W/cm**3" % ( j_total * 1.0e-6 )
        j_total = rsm.radiative_integrated_emission_for_gas_state( Q, False )
        print "Spectrally unresolved integrated emission: j_total = %0.2f W/cm**3" % ( j_total * 1.0e-6 )
        del X
        
        print "Testing read from file..."
        X = CoeffSpectra()
        X.read_from_file("coefficient_spectra.txt",0,-1);
        j_total = X.write_to_file( "coefficient_spectra.txt", WAVELENGTH );
        
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
        print "Total radiative flux from a %f cm slab: q_total = %0.2f W/cm**2 (%e W/m**2)" % ( (s1-s0)*100, q_total*1.0e-4, q_total )
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

        print  "testing quick_solve_for_divq_with_binning(OPACITY_BINNING, N_bins=10)"
        q_rad = TS.quick_solve_for_divq_with_binning( OPACITY_BINNING, 10 )
        print  "q_rad = %e W/m2 " % q_rad
        print  "testing quick_solve_for_divq_with_binning(OPACITY_BINNING, N_bins=100)"
        q_rad = TS.quick_solve_for_divq_with_binning( OPACITY_BINNING, 100 )
        print  "q_rad = %e W/m2 " % q_rad
        print  "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=10)"
        q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 10 )
        print  "q_rad = %e W/m2 " % q_rad
        print  "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=95)"
        q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 95 )
        print  "q_rad = %e W/m2 " % q_rad
        print  "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=950)"
        q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 950 )
        print  "q_rad = %e W/m2 " % q_rad
        print  "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=9500)"
        q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 9500 )
        print  "q_rad = %e W/m2 " % q_rad
        print  "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=95000)"
        q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 95000 )
        print  "q_rad = %e W/m2 " % q_rad

        del TS
        
        print  "testing quick_solve_for_divq_in_blocks(OPACITY_BINNING, N_blocks=1)"
        TS = TS_data(rsm,nslabs,T0,T1,True)
        for islab in range(nslabs):
            divq = new_doublep()
            TS.set_rad_point(islab,Q,divq,s0+0.5*ds+islab*ds,ds)
        q_rad = TS.quick_solve_for_divq_in_blocks( 1 )
        print  "q_rad = %e W/m2 " % q_rad
        print  "testing quick_solve_for_divq_in_blocks(N_blocks=10)"
        q_rad = TS.quick_solve_for_divq_in_blocks( 10 )
        print  "q_rad = %e W/m2 " % q_rad
        
        del TS

    del gm, rsm, Q
    
    print "done."
    
if __name__ == '__main__':
    main()
