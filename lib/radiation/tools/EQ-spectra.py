#! /usr/bin/python

import sys
from gaspy import *
from radpy import *
from cfpylib.gasdyn.cea2_gas import *
from cfpylib.util.YvX import *
from time import time
from getopt import getopt, GetoptError

# sampling scale factor
f_s = 1.0 

from gaspy import *

longOptions = ["help", "input-file="]

def printUsage():
    print ""
    print "Usage: eq-spectra.py [--help] [--input-file=<fileName>]"
    print "e.g. eq-spectra.py --input-file='eqs.py'"
    print ""
    return
    
class EqSpectraInputData(object):
    """Python class to organize the global data for the EQ-spectra calculation.

    The user's script does not create this object but rather just alters
    attributes of the global object.
    """
    __slots__ = ['rad_model_file', 'gas_model_file', 'species_list', 'mole_fractions', 
                 'mass_fractions', 'shock_speed', 'gas_pressure', 'gas_temperature', 
                 'tube_width', 'apparatus_fn', 'Gaussian_HWHM', 'Lorentzian_HWHM', 
                 'sampling_rate', 'problem', 'planck_spectrum', 'show_plots' ]
    def __init__(self):
        self.rad_model_file = "rad-model.lua"
        self.gas_model_file = "gas-model.lua"
        self.species_list = [ 'Ar', 'Ar_plus', 'e_minus' ]
        self.mole_fractions = None
        self.mass_fractions = None
        self.shock_speed = 0
        self.gas_pressure = 0
        self.gas_temperature = 0
        self.tube_width = 0
        self.apparatus_fn = None
        self.Gaussian_HWHM = 0
        self.Lorentzian_HWHM = 0
        self.sampling_rate = 1
        self.problem = "shock"
        self.planck_spectrum = False
        self.show_plots = True
    
def parseInputData(input_data):
    # parse the input data object
          
    print "Setting up a radiation model from input file", input_data.rad_model_file 
    rsm = create_radiation_spectral_model( input_data.rad_model_file )
    
    species = input_data.species_list
    # make the gas-model
    create_gas_file( "thermally perfect gas", species, input_data.gas_model_file )
    print "Setting up a gas model from input file", input_data.gas_model_file 
    gm = create_gas_model(input_data.gas_model_file)
    nsp = gm.get_number_of_species()
    ntm = gm.get_number_of_modes()
    
    # species mass-fractions or mole-fractions
    if input_data.mass_fractions!=None:
        massf_inf = input_data.mass_fractions
    else:
        massf_inf = convert_molef2massf(input_data.mole_fractions,gm.M())
    # print "massf_inf = ", massf_inf        

    # shock speed
    Us = input_data.shock_speed
    
    # gas pressure
    p_inf = input_data.gas_pressure

    # gas temperature
    T_inf = input_data.gas_temperature
    
    # tube width
    tube_D = input_data.tube_width

    # name of the desired apparatus function model
    apparatus_fn = input_data.apparatus_fn
    
    # Gaussian half-width half maximum of the spectrometer apparatus function
    gamma_G = input_data.Gaussian_HWHM
            
    # Lorentzian half-width half maximum of the spectrometer apparatus function
    gamma_L = input_data.Lorentzian_HWHM    
            
    # sampling output for the plotting-program-readable output
    nu_sample = input_data.sampling_rate
    
    # the problem type
    problem = input_data.problem
    
    # planck spectrum request
    planck_spectrum = input_data.planck_spectrum
    
    show_plots = input_data.show_plots
    
    return rsm, gm, species, nsp, ntm, massf_inf, Us, p_inf, T_inf, tube_D, apparatus_fn, gamma_G, gamma_L, nu_sample, problem, planck_spectrum, show_plots

def main():
    #
    try:
        userOptions = getopt(sys.argv[1:], [], longOptions)
    except GetoptError, e:
        print "One (or more) of your command-line options was no good."
        print "    ", e
        printUsage()
        sys.exit(1)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(0)
    #
    input_file = uoDict.get("--input-file", "none")
    #
    # create the input data instance
    input_data = EqSpectraInputData()
    # 
    # parse the input options
    execfile(input_file)
    rsm, gm, species, nsp, ntm, massf_inf, Us, p_inf, T_inf, tube_D, apparatus_fn, gamma_G, gamma_L, nu_sample, problem, planck_spectrum, show_plots = parseInputData(input_data)

    # setup the reactants list
    reactants = make_reactants_dictionary( species )
    for sp in species:
        _sp = sp.replace("_plus","+").replace("_minus","-")
        reactants[_sp] = massf_inf[gm.get_isp_from_species_name(sp)]
    
    massf_sum = 0.0
    for massf in reactants.values():
        massf_sum += massf
    for sp in reactants.keys():
        reactants[sp] /= massf_sum 
    print reactants

    # solve the post-shock equilibrium radiation problem
    Q = Gas_data(gm)
    Q.p = p_inf
    for itm in range(ntm): Q.T[itm] = T_inf
    for isp,sp in enumerate(species): Q.massf[isp] = get_species_composition(sp,reactants)
    # firstly without the onlyList:
    cea = Gas( reactants, with_ions=get_with_ions_flag(species), trace=0 )
    cea.set_pT(p_inf,T_inf,transProps=False)
    if problem=="shock":
        cea.shock_process( Us )
    print "unfiltered species composition: ", cea.species
    del cea 
    cea = Gas( reactants, onlyList=reactants.keys(), with_ions=get_with_ions_flag(species), trace=1.0e-20 )
    cea.set_pT(p_inf,T_inf,transProps=False)
    if problem=="shock":
        cea.shock_process( Us )
    print "filtered species composition: ", cea.species
    #over-write provided initial mass-fractions
    Q.rho = cea.rho
    for itm in range(ntm): Q.T[itm] = cea.T
    for isp,sp in enumerate(species):
        Q.massf[isp] = get_species_composition(sp,cea.species)
            
    gm.eval_thermo_state_rhoT(Q)
    print "computed equlibrium state with filtered species: "
    Q.print_values(False)

    if "e_minus" in species:
        N_elecs = ( Q.massf[gm.get_isp_from_species_name("e_minus")]*Q.rho/RC_m_SI*1.0e-6 )
        N_total = ( ( Q.p - Q.p_e ) / RC_R_u / Q.T[0] + Q.p_e / RC_R_u / Q.T[-1] ) * RC_Na * 1.0e-6
        print "electron number density = %e cm-3" % ( N_elecs )
        print "total number density = %e cm-3" % N_total
        print "ionization fraction = %e" % ( N_elecs / N_total )

    # perform LOS calculation
    LOS = LOS_data( rsm, 1 )
    Q_rE_rad = new_doublep()
    t0 = time()
    LOS.set_rad_point(0,Q,Q_rE_rad,tube_D*0.5,tube_D)
    t1 = time()
    print "Wall time = %f seconds" % ( t1-t0 ) 
    LOS.write_point_to_file(0,"coefficient_spectra.txt")
    S = SpectralIntensity( rsm )
    I_total = LOS.integrate_LOS( S )
    print "I_total = %0.3e W/m2" % I_total
    # initialise apparatus function
    if apparatus_fn=="Voigt":
        A = Voigt(gamma_L, gamma_G, nu_sample)
    elif apparatus_fn=="SQRT_Voigt":
        A = SQRT_Voigt(gamma_L, gamma_G, nu_sample)
    elif apparatus_fn=="none" or apparatus_fn=="None" or apparatus_fn==None:
        A = None
    else:
        print "Apparatus function with name: %s not recognised." % apparatus_fn
        sys.exit()
    
    # apply apparatus function
    if A!=None:
        S.apply_apparatus_function(A)
    S.write_to_file("intensity_spectra.txt" ) 
    
    if show_plots:
        IvW = YvX("intensity_spectra.txt" )
        IvW.plot_data(xlabel="Wavelength, lambda (nm)", ylabel="Intensity, I (W/m2-m-sr)", new_plot=True, show_plot=True, include_integral=False, logscale_y=True )
        del IvW
    
    # check if the planck intensity spectrum has been requested
    if planck_spectrum:
        del S
        S = SpectralIntensity(rsm,cea.T)
        S.write_to_file("planck_intensity_spectra.txt" )
    
    del S
    
    del rsm, gm
    
    print "done."
    
if __name__ == '__main__':
    main()
