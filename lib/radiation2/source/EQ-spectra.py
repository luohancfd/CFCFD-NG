#! /usr/bin/python

import sys
from gaspy import *
from librad2 import *
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
    print "e.g. eq-spectra.py --input-file='eqs.inp'"
    print ""
    return
    
def parseInputFile(input_file):
    # parse the input file
    ifile = open(input_file,"r")
    lines = ifile.readlines()
    ifile.close()
    if len(lines)!=11:
        print "Input file %s does not follow the expected format: " % input_file
        print "rad-model-file: <radFile>"
        print "   species-list: <speciesList>"
        print " mole-fractions: <moleFractionsList>"
        print "     --- OR ---"
        print " mass-fractions: <massFractionsList>"
        print "    shock-speed: <shockSpeed> <units>"
        print "   gas-pressure: <gasPressure> <units>"
        print "gas-temperature: <gasTemperature> <units>"
        print "     tube-width: <tubeWidth> <units>"
        print "   Apparatus-fn: <apparatus_fn>"
        print "  Gaussian-HWHM: <gamma_G> <units>"
        print "Lorentzian-HWHM: <gamma_L> <units>"
        print "  sampling-rate: <nu_sample>"
        sys.exit()
          
    # first line is the radiation mode filename
    tks = lines[0].split()
    rmodel_file = tks[1]
    # make the radiation-model
    print "Setting up a radiation model from input file", rmodel_file 
    rsm = create_radiation_spectral_model( rmodel_file )
    
    # second line is the species list
    tks = lines[1].split()
    species = tks[1:]
    # make the gas-model
    gmodel_file = "gas-model.lua"
    create_gas_file( "thermally perfect gas", species, gmodel_file )
    print "Setting up a gas model from input file", gmodel_file 
    gm = create_gas_model("gas-model.lua")
    nsp = gm.get_number_of_species()
    ntm = gm.get_number_of_modes()
    
    # third line are the species mass-fractions or mole-fractions
    tks = lines[2].split()
    if tks[0]=="mass-fractions:":
        massf_inf = []
        for tk in tks[1:]:
            massf_inf = float(tk)
    elif tks[0]=="mole-fractions:":
        molef_inf = []
        for tk in tks[1:]:
            molef_inf.append( float(tk) )
        massf_inf = convert_molef2massf(molef_inf,gm.M())
    # print "massf_inf = ", massf_inf        

    # forth line is the shock speed
    tks = lines[3].split()
    Us = float(tks[1])
    if len(tks)==3:
        if tks[2]=="km/s": Us *= 1000.
        elif tks[2]=="m/s": Us *= 1.0
        else:
            print "Shock speed units: %s not understood" % tks[2]
            sys.exit()
    
    # fifth line is the gas pressure
    tks = lines[4].split()
    p_inf = float(tks[1])
    if len(tks)==3:
        if tks[2]=="Torr" or tks[2]=="torr": p_inf *= 133.3333
        elif tks[2]=="Pa" or tks[2]=="pa": p_inf *= 1.0
        else:
            print "Pressure units: %s not understood" % tks[2]
            sys.exit()
    
    # sixth line is the gas temperature
    tks = lines[5].split()
    T_inf = float(tks[1])
    
    # seventh line is the tube width
    tks = lines[6].split()
    tube_D = float(tks[1])
    if len(tks)==3:
        if tks[2]=="cm": tube_D /= 100.0
        elif tks[2]=="mm": tube_D /= 1000.0
        elif tks[2]=="m": tube_D *= 1.0
        else:
            print "Tube width units: %s not understood" % tks[2]
            sys.exit()

    # eighth line is the name of the desired apparatus function model
    tks = lines[7].split()
    apparatus_fn = tks[1]
    
    # ninth line is the Gaussian half-width half maximum of the spectrometer apparatus function
    tks = lines[8].split()
    gamma_G = float(tks[1])
    if len(tks)==3:
        if tks[2]=="Ang" or tks[2]=="ang": gamma_G *= 1.0
        elif tks[2]=="nm": gamma_G *= 10.0
        else:
            print "FWHM units: %s not understood" % tks[2]
            sys.exit()
            
    # tenth line is the Lorentzian half-width half maximum of the spectrometer apparatus function
    tks = lines[9].split()
    gamma_L = float(tks[1])
    if len(tks)==3:
        if tks[2]=="Ang" or tks[2]=="ang": gamma_L *= 1.0
        elif tks[2]=="nm": gamma_L *= 10.0
        else:
            print "FWHM units: %s not understood" % tks[2]
            sys.exit()     
            
    # eleventh line is the sampling output for the plotting-program-readable output
    tks = lines[10].split()
    nu_sample = int(tks[1])  
    
    return rsm, gm, species, nsp, ntm, massf_inf, Us, p_inf, T_inf, tube_D, apparatus_fn, gamma_G, gamma_L, nu_sample

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
    # parse the input option
    rsm, gm, species, nsp, ntm, massf_inf, Us, p_inf, T_inf, tube_D, apparatus_fn, gamma_G, gamma_L, nu_sample = parseInputFile(input_file)

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
    cea.shock_process( Us )
    # print "unfiltered species composition: ", cea.species
    del cea 
    cea = Gas( reactants, onlyList=reactants.keys(), with_ions=get_with_ions_flag(species), trace=1.0e-20 )
    cea.set_pT(p_inf,T_inf*10,transProps=False)
    cea.T = T_inf
    cea.shock_process( Us )
    # print "filtered species composition: ", cea.species
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
    print "I_total = ", I_total
    # initialise apparatus function
    if apparatus_fn=="Voigt":
        A = Voigt(gamma_L, gamma_G, nu_sample)
    elif apparatus_fn=="SQRT_Voigt":
        A = SQRT_Voigt(gamma_L, gamma_G, nu_sample)
    elif apparatus_fn=="none" or apparatus_fn=="None":
        A = None
    else:
        print "Apparatus function with name: %s not recognised." % apparatus_fn
        sys.exiit()
    
    # apply apparatus function
    if A!=None:
        S.apply_apparatus_function(A)
    S.write_to_file("intensity_spectra.txt" ) 
    
    IvW = YvX("intensity_spectra.txt" )
    IvW.plot_data(xlabel="Wavelength, lambda (nm)", ylabel="Intensity, I (W/m2-m-sr)", new_plot=True, show_plot=True, include_integral=False, logscale_y=True )
    
    del IvW
    
    del rsm, gm
    
    print "done."
    
if __name__ == '__main__':
    main()
