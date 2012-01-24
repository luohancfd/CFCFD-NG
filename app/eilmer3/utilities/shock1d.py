#! /usr/bin/env python
## \file shock1d.py
## \ingroup basic_gas_dyn
## \brief Python program to compute the flow properties across a 1D shock.
##
## \author P.Jacobs
## \version 09-May-2007 adapted from shock1d.c
## \version 18-Jan-2011 ported to use Rowan's latest gas models.

import sys, os
from getopt import getopt
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from libprep3 import *
from cfpylib.gasdyn.ideal_gas_flow import m2_shock, r2_r1, p2_p1, u2_u1, T2_T1

shortOptions = ""
longOptions = ["help", "demo", "gas=", "file=", "species=", "p=", "T=", "massf=", "us="]

def printUsage():
    print ""
    print "Usage: shock1d.py" + \
          " [--help]" + \
          " [--demo]" + \
          " [--gas=<gas-model-name>]" + \
          " [--file=<gas-model-file-name>]" + \
          " [--species=<list of species names>]" + \
          " [--p=<pressure,Pa>]" + \
          " [--T=<temperature,K>]" + \
          " [--massf=<list of mass-fractions>]" + \
          " [--us=<shock-speed,m/s>]"
    print "Note that all of the arguments are optional."
    print "Mass-fractions should be specified as a comma-separated list with no spaces."
    print ""
    return

def select_gas_model(model=None, species=None, fname=None):
    """
    Selects a gas model for the calculation.

    Input:
    model   : (string) name of the gas model as shown in the list below.
    species : list of species names (strings).
    fname   : (string) name of the gas-model file.

    The new gas models are configured by stand-alone files.
    If you already have a gas-model.lua file already set up,
    give its name af fname.
    If you have not already set up the gas-model.lua file,
    this function is provided as a simple but limited means to do so.

    Look-up-table (LUT) gas models and LUT_plus_composite cannot be
    created directly with this function.
    If you want a single LUT gas model, just set up the LUT table 
    file externally and supply the name of that file as fname.
    If you want a LUT-plus-composite gas model, set up the LUT table
    externally and then set up the rest of the composite gas model
    using create_gas_file(), which has the capability of prepending 
    the LUT gas species.
    """
    if fname is None or fname is "":
        # Help the user to set up the gas-model file.
        fname = "gas-model.lua"
        if model == None:
            print "select_gas_model():"
            print "    A gas 'model' or 'fname' must be specified."
            print "    Bailing out!"
            sys.exit(1)
        if species == None:
            print "select_gas_model():"
            print "    When setting up a gas model, a list of species must be specified."
            print "    Bailing out!"
            sys.exit(1)
        create_gas_file(model, species, fname)
    gmodel = set_gas_model_ptr(create_gas_model(fname))
    nsp = gmodel.get_number_of_species()
    return [ gmodel.species_name(isp) for isp in range(nsp) ]


if __name__ == '__main__':
    print "shock1d: Compute properties across a 1D shock."
    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(0)
    #
    if uoDict.has_key("--demo"):
        print "Begin demonstration using a look-up table for air."
        print "Example from Figure 14.2 in Anderson's hypersonics text."
        gas_model_file = "cea-lut-air.lua.gz"
        if not os.access(gas_model_file, os.R_OK):
            print "Prepare CEA look-up table."
            os.system("build-cea-lut --case=air")
        gas_model_name = None  # all info is in the .lua file
        gas_species = None     # ditto
        massf = [1.0,]         # mass fractions
        p1 = 10.1e3            # initial pressure, Pa
        T1 = 225.0             # initial temperature, degrees K
        us = 6.0e3             # shock speed, m/s
    else:
        # Set up the problem from command-line arguments
        # The default gas model is ideal air.
        gas_model_name = uoDict.get("--gas", "ideal gas")
        gas_model_file = uoDict.get("--file", "")
        if uoDict.has_key("--species"):
            species_string = uoDict.get("--species", "air")
            gas_species = [item.strip() for item in species_string.split(',')]
        else:
            gas_species = ["air",]
        p1 = float(uoDict.get("--p", "100000.0"))
        T1 = float(uoDict.get("--T", "296.0"))
        if uoDict.has_key("--massf"):
            massf_string = uoDict.get("--massf", "1.0")
            massf = [float(item) for item in massf_string.split(',')]
            if abs(sum(massf) - 1.0) > 1.0e-6:
                print "Mass fractions sum to", sum(massf), "but should sum to 1.0"
        else:
            massf = [1.0,]
        us = float(uoDict.get("--us", "3000.0"))
    #    
    species_list = select_gas_model(model=gas_model_name, 
                                    species=gas_species,
                                    fname=gas_model_file)
    print "species_list=", species_list
    gmodel = get_gas_model_ptr()
    nsp = gmodel.get_number_of_species()
    nmodes = gmodel.get_number_of_modes()
    if nmodes > 1:
        print "This program uses only mode0 energy."
        print "If this is a problem for your gas model,"
        print "some more code needs to be written."
        sys.exit(1)

    q1 = Gas_data(gmodel)
    q1.p = p1
    for i in range(nmodes): q1.T[i] = T1
    for i in range(nsp): q1.massf[i] = massf[i]

    gmodel.eval_thermo_state_pT(q1)
    gama = gmodel.gamma(q1)
    R = gmodel.R(q1)
    print "Preshock gas state:"
    print "    p=", q1.p, " T=", q1.T[0]
    print "    rho=", q1.rho, " e=", q1.e[0], " a=", q1.a
    print "    mass-fractions=", massf
    print "    gas-gamma=", gama, " R=", R
    #
    M1 = us / q1.a
    print "Shock-speed=", us, " Mach-number=", M1
    #
    print "Post-shock condition, assuming ideal gas:"
    M2 = m2_shock(M1, gama)
    p2 = p2_p1(M1, gama) * q1.p
    rho2 = r2_r1(M1, gama) * q1.rho
    T2 = T2_T1(M1, gama) * q1.T[0]
    u2 = u2_u1(M1, gama) * us # gas speed in shock frame
    ug = us - u2  # gas speed in lab frame
    print "    M2=", M2, " u2=", u2, " ug=", ug
    print "    p2=", p2, " T2=", T2, " rho2=", rho2
    #
    print "Post-shock condition, using full gas model:"
    q2 = Gas_data(gmodel)
    q2.p = p2
    for i in range(nmodes): q2.T[i] = T2
    for i in range(nsp): q2.massf[i] = massf[i]
    gmodel.eval_thermo_state_pT(q2)
    #
    q3 = Gas_data(gmodel)
    q3.p = p2
    for i in range(nmodes): q3.T[i] = T2
    for i in range(nsp): q3.massf[i] = massf[i]
    gmodel.eval_thermo_state_pT(q3)

    e1 = q1.e[0]
    rho1 = q1.rho
    V1 = us
    temp1 = p1 + rho1 * V1 * V1  # momentum
    H1 = e1 + p1 / rho1 + 0.5 * V1 * V1  # total enthalpy
    rho_delta = 1.0
    e_delta = 1.0
    tol = 1.0e-8
    # Update the estimates using the Newton-Raphson method.
    for count in range(20):
        rho_save = q2.rho
        e_save = q2.e[0]
        p_save = q2.p
        #
        p2 = p_save
        e2 = e_save
        rho2 = rho_save
        r2r1 = rho2 / rho1
        f1 = temp1 - p2 - rho1 * rho1 * V1 * V1 / rho2
        f2 = H1 - e2 - p2 / rho2 - 0.5 * V1 * V1 / (r2r1 * r2r1)
        f1_save = f1
        f2_save = f2
        #
        # Use finite differences to compute the Jacobian
        d_rho = rho_save * 0.01
        d_e = e_save * 0.01
        #
        rho2 = rho_save + d_rho
        e2 = e_save
        q3.rho = rho2
        q3.e[0] = e2
        gmodel.eval_thermo_state_rhoe(q3)
        p2 = q3.p
        f1 = temp1 - p2 - rho1 * rho1 * V1 * V1 / rho2
        f2 = H1 - e2 - p2 / rho2 - 0.5 * V1 * V1 / (r2r1 * r2r1)
        A = (f1 - f1_save) / d_rho
        C = (f2 - f2_save) / d_rho
        #
        rho2 = rho_save
        e2 = e_save + d_e
        q3.rho = rho2
        q3.e[0] = e2
        gmodel.eval_thermo_state_rhoe(q3)
        p2 = q3.p
        f1 = temp1 - p2 - rho1 * rho1 * V1 * V1 / rho2
        f2 = H1 - e2 - p2 / rho2 - 0.5 * V1 * V1 / (r2r1 * r2r1)
        B = (f1 - f1_save) / d_e
        D = (f2 - f2_save) / d_e
        #
        # Invert Jacobian and multiply.
        det = A * D - B * C
        rho_delta = (D * f1_save - B * f2_save) / det
        e_delta = (-C * f1_save + A * f2_save) / det
        #
        rho_new = rho_save - rho_delta
        e_new   = e_save - e_delta
        #
        q2.rho = rho_new
        q2.e[0] = e_new
        gmodel.eval_thermo_state_rhoe(q2)
        #
        # Check convergence.
        if abs(rho_delta) < tol and abs(e_delta) < tol: break
    print "    count=", count, " rho_delta=", rho_delta, " e_delta=", e_delta
    #
    # Back-out velocities via continuity.
    V2 = V1 * q1.rho / q2.rho
    M2 = V2 / q2.a
    ug = V1 - V2
    print "    M2=", M2, " u2=", V2, " ug=", ug, "Mg=", ug/q2.a
    print "    p2=", q2.p, " T2=", q2.T[0]
    print "    rho2=", q2.rho, " e2=", q2.e[0], " a2=", q2.a

    if uoDict.has_key("--demo"):
        print "Expected T2/T1= 35 (approximately)"
        print "Computed T2/T1=", q2.T[0]/q1.T[0]

print "Done."
