#! /usr/bin/python

from gaspy import *
from radpy import *
from cfpylib.gasdyn.cea2_gas import *
from numpy import *
import sys
   
def main():
    print "--------------------------------------------------------------------------"
    print "|       RADIATION OF HIGH TEMPERATURE GASES TESTCASE TC4 LEVEL 1         |"
    print "--------------------------------------------------------------------------"
    print "*   - Equilibrium N2-CH4 plasma torch emission                           *"
    print "*   - Optically thick line-of-sight with intensity attenuation           *"
    print "*   - Radiating species considered are  C2, CN, N2                       *"
    print "*   - Implementing photaura radiation model from librad                  *"
    print "--------------------------------------------------------------------------\n"

    # 1. Setup the gas model
    print "Creating the gas-model for a N2-CH4 gas mix"
    gm = create_gas_model("gas-model.lua")
    nsp = gm.get_number_of_species()
    ntm = gm.get_number_of_modes()    
    # make the species list
    species = []
    for isp in range(nsp):
        species.append( gm.species_name(isp) )
    # initialise a new gas-data instance
    Q = Gas_data(gm)
    
    # 2. Setup the radiation model
    print "Setting up a photaura spectral radiation model"
    rsm = create_radiation_spectral_model("TC4-radiators.lua")
    
    # 3. Setup calculation control parameters
    HWHM = 11.0 / 2
    print "Using a spectrometer HWHM of %0.2f Angstroms" % HWHM
    nu_skip = 100
    do_LOS_integration = True
    extract_emissivities = True
    n_slabs = 5 # there are approximately 110 temperature data points, so 110 should be the maximum value used here
    x_list= [ 1.24e-3, 3.82e-3, 8.43e-3 ] # list of locations to extract the emissivities at

    # 4. Setup the CEA equilibrium calculator 
    p_plasma = 30.0e3     # Pa
    Q.p = p_plasma
    CH4_mw = gm.molecular_weight( gm.get_isp_from_species_name( "CH4" ) )
    N2_mw = gm.molecular_weight( gm.get_isp_from_species_name( "N2" ) )
    # 98.1% N2 by volume
    CH4_mf = 0.019 * CH4_mw / ( 0.019 * CH4_mw + 0.981 * N2_mw )
    N2_mf  = 0.981 * N2_mw / ( 0.019 * CH4_mw + 0.981 * N2_mw )
    #
    reactants = make_reactants_dictionary( species )
    reactants["CH4"] = CH4_mf
    reactants["N2"] = N2_mf
    cea = Gas( reactants, onlyList=reactants.keys(), with_ions=True, trace=1.0e-10 )
    
    # 5. Perform line-of-sight integration
    if do_LOS_integration:
        # Initialise line-of-sight class instance (f_res and radial limit previously defined)
        T_i = 0.0; T_f = 0.0
        Plasma = LOS_data(rsm, n_slabs, T_i, T_f)

        # provided temperature profile/s:
        profile_type = "nominal"
        if profile_type == "nominal":
            x_points = [ -10.64e-3, -10.44e-3, -10.22e-3, -10.05e-3, -9.8e-3, -9.64e-3, -9.4e-3, -9.23e-3, -9.03e-3, -8.81e-3, -8.651e-3, -8.436e-3, -8.211e-3, -8.011e-3, -7.836e-3, -7.611e-3, -7.426e-3, -7.221e-3, -7.021e-3, -6.841e-3, -6.616e-3, -6.441e-3, -6.226e-3, -6.001e-3, -5.821e-3, -5.651e-3, -5.421e-3, -5.236e-3, -5.016e-3, -4.839e-3, -4.601e-3, -4.406e-3, -4.216e-3, -4.016e-3, -3.821e-3, -3.621e-3, -3.426e-3, -3.221e-3, -2.996e-3, -2.831e-3, -2.631e-3, -2.415e-3, -2.236e-3, -2.031e-3, -1.826e-3, -1.621e-3, -1.426e-3, -1.241e-3, -1.021e-3, -0.851e-3, -0.605e-3, -0.431e-3, -0.353e-3, -0.241e-3, -0.163e-3, -0.045e-3, 0.0e-3,  0.045e-3, 0.163e-3, 0.241e-3, 0.353e-3, 0.431e-3, 0.605e-3, 0.851e-3, 1.021e-3, 1.241e-3, 1.426e-3, 1.621e-3, 1.826e-3, 2.031e-3, 2.236e-3, 2.415e-3, 2.631e-3, 2.831e-3, 2.996e-3, 3.221e-3, 3.426e-3, 3.621e-3, 3.821e-3, 4.016e-3, 4.216e-3, 4.406e-3, 4.601e-3, 4.839e-3, 5.016e-3, 5.236e-3, 5.421e-3, 5.651e-3, 5.821e-3, 6.001e-3, 6.226e-3, 6.441e-3, 6.616e-3, 6.841e-3, 7.021e-3, 7.221e-3, 7.426e-3, 7.611e-3, 7.836e-3, 8.011e-3, 8.211e-3, 8.436e-3, 8.651e-3, 8.810e-3, 9.030e-3, 9.230e-3, 9.400e-3, 9.640e-3, 9.800e-3, 10.05e-3, 10.22e-3, 10.44e-3, 10.64e-3]
            T_points = [  3041.0,    3093.0,   3150.0,     3194.0,    3258.0,  3300.0,   3362.0,  3406.0,   3458.0,   3515.0,   3563.0,    3634.0,     3699.0,    3755.0,    3816.0,    3863.0,    3912.0,    3959.0,    4002.0,    4058.0,   4104.0,     4164.0,    4231.0,    4286.0,   4337.0,     4401.0,    4449.0,    4499.0,   4535.0,    4576.0,     4605.0,    4629.0,    4653.0,    4678.0,    4706.0,    4735.0,    4768.0,    4806.0,    4834.0,    4867.0,    4899.0,    4923.0,    4946.0,    4961.0,    4967.0,   4968.0,    4965.0,     4959.0,    4953.0,   4946.0,    4945.0,     4940.0,    4940.0,    4937.0,   4937.0,    4936.0,   4936.0,  4936.0,   4937.0,   4937.0,    4940.0,   4940.0,  4945.0,    4946.0,   4953.0,    4959.0, 4965.0,    4968.0,   4967.0,   4961.0,   4946.0,   4923.0,   4899.0,   4867.0,    4834.0, 4806.0,     4768.0,   4735.0, 4706.0,    4678.0,    4653.0,  4629.0,   4605.0,    4576.0,  4535.0,    4499.0,  4449.0,   4401.0,   4337.0,   4286.0,   4231.0,   4164.0,   4104.0,   4058.0,   4002.0,   3959.0,   3912.0,   3863.0,   3816.0,   3755.0,    3699.0, 3634.0,    3563.0,   3515.0,   3458.0,   3406.0,   3362.0,    3300.0, 3258.0,    3194.0,    3150.0,   3093.0,  3041.0 ]
        else:
            print "profile_type %s is not defined, exiting program..." % profile_type
            sys.exit()
    
        # ensure we have an odd number of slabs
        if n_slabs%2 != 1:
            print "Symetrical geometry, need an odd number of slabs"
            sys.exit()

        # discretise into n_slabs (defined at top of page)
        def slabs( x ):
            return interp( x, x_points, T_points )
        # slabs.plot_spline(show_plot=True, include_integral=False)
        slab_width = ( x_points[-1] - x_points[0] )/n_slabs
    
        # write the temperature profile to file
        if 0:
            ofile = open( "spline.txt", 'w' )
            dx = 1.0e-4
            x = x_points[0]
            while x<=x_points[-1]:
                ofile.write("%e \t %e\n" % ( x, slabs(x) ) )
                x += dx
            ofile.close()
    
        # calculate and store the spectral coefficients through the plume
        divqs = []
        for islb in range(n_slabs):
            x = (islb + 0.5) * slab_width + x_points[0]
            # find what temp fits this x through linear interpolation
            T = slabs( x )
            # create the divq pointer
            divqs.append( new_doublep() )
            # get on with preparing the gas-state
            print "Preparing slab %d of %d: x = %0.2f mm, T = %0.1f K" % ( islb+1, n_slabs, x*1.e3, T )
            # only do unique calculations up to centre slab
            if islb<(n_slabs+1)/2:
                cea.set_pT(p_plasma,T)
                for itm in range(ntm): Q.T[itm] = T
		for isp,sp in enumerate(species):
		    Q.massf[isp] = get_species_composition(sp,cea.species)
                gm.eval_thermo_state_pT(Q)
                # Q.print_values(False)
                Plasma.set_rad_point(islb,Q,divqs[islb],x,slab_width)
            else:
                # partner = centre_slab        -  distance from centre slab
                ipslb = ( (n_slabs+1)/2 - 1 ) - ( islb - ( (n_slabs+1)/2 - 1 ) )
                print "islb %d has partner slab %d" % ( islb + 1, ipslb + 1 )
                Plasma.clone_rad_point(ipslb,islb,divqs[islb],x,slab_width)
                print "nu.size() = ", Plasma.get_rpoint_pointer(islb).X_.nu.size()
    
        # create SpectralIntensity instance to hold the spectra
        S = SpectralIntensity( rsm )
	
        print "\nperforming spatial integration through plasma torch..."
        I_total = Plasma.integrate_LOS(S)
    
        # make a copy to retain unsmeared spectrum
        S_copy = SpectralIntensity(S)
        # create the Gaussian apparatus function
        A = Voigt(0.0, HWHM, nu_skip)
        # apply apparatus function
        S.apply_apparatus_function( A )
    
        # write to file
        outfile = "plasma-transported-intensity.txt"
        I_total = S.write_to_file( outfile )
        print "I_total after plasma transport = %0.3e W/m**2-sr" % ( I_total )
        
        # compute 'optically thin' intensity
        for inu in range(S_copy.nu.size()):
            S_copy.I_nu[inu] = 0.0
        for islb in range(n_slabs):
            for inu in range(S_copy.nu.size()):
                S_copy.I_nu[inu] += Plasma.get_rpoint_pointer(islb).X_.j_nu[inu] * slab_width
        # apply apparatus function
        S_copy.apply_apparatus_function( A )
        # write to file
        outfile = "plasma-transported-intensity-optically-thin.txt"
        I_total = S_copy.write_to_file( outfile )
        print "I_total after plasma transport (optically thin) = %0.3e W/m**2-sr" % ( I_total )
        
        # clean up some c++ data structures
        del S, S_copy
        for divq in divqs:
            delete_doublep( divq )

    # 6. extract emissivities for comparison with Abel inversions    
    if extract_emissivities:
        for x in x_list:
            # gas state at this location
            T = slabs(x)
            cea.set_pT(p_plasma,T)
            for itm in range(ntm): Q.T[itm] = T
            for isp,sp in enumerate(species):
                Q.massf[isp] = get_species_composition(sp,cea.species)
            gm.eval_thermo_state_pT(Q)
            
            # check the CN and C2 mole-fractions
            M = vectord(nsp)
            molef = vectord(nsp)
            for isp in range(nsp):
                M[isp] = gm.molecular_weight(isp)
            convert_massf2molef(Q.massf,M,molef)
            print "x = %e m" % x
            print "X[CN] = ", molef[species.index("CN")]
            print "X[C2] = ", molef[species.index("C2")]
            
            # calculate emission and absorption coefficient spectra
            X = CoeffSpectra(rsm)
            rsm.radiative_spectra_for_gas_state(Q,X)
            
            # apply apparatus function to the coefficient spectra
            X.apply_apparatus_function( A )
            X.integrate_emission_spectra()
            
            # write to file
            X.write_to_file("emissivities_at_%d_microns.txt" % int(x*1000000) )
            del X,S

    # delete the gas and spectral radiation models
    del gm, rsm

    print "\nTC4 testcase validation complete."
    
if __name__ == '__main__':
    main()
