#! /usr/bin/env python

## \file    reactor.py
## \author  Brendan O'Flaherty
## 
## \version 24 July 2008 -- File created
## 
## \brief   A reactor that uses coupled operations.
## ---------------------------------------------------------

import os
import sys
import ConfigParser
from getopt import getopt
from time import localtime

from numpy import array, arange, zeros

sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *
from scipy.integrate import *

# ----------------------------------------------------------

def const_v_fixed_mass(y, t, Q, g, r, dt_chem, testFlag):
    # Y - mass fractions (kg/kg)
    # w - concentration (mol/m**3)
    # dwdt - molar production rate (mol/s)
    # M - molecular weights (kg/mol)
    
    nsp = g.get_number_of_species()

    # unpack array
    Q.T[0] = y[0] # temperature
    Q.p = y[1] # pressure
    Y = y[2:nsp+2] # mass fraction

    dYdt = zeros(nsp)

    dwdt = vectord(nsp*[0.0])
    w = vectord(nsp*[0.0])
    
    # update gas
    for isp in range(nsp):
        Q.massf[isp] = Y[isp]

    g.eval_thermo_state_pT(Q) # concs _are_ updated here

    # copy concs to w
    w = convert_massf2conc(Q.rho, Q.massf, g.M())

    # get rates-of-change
    dwdt = r.rate_of_change_py(Q)

    sum_ew = 0.0
    sum_wcv = 0.0
    for i in range(nsp):
        e_i = g.internal_energy(Q, i)*g.M()[i]
        sum_ew += e_i*dwdt[i]
        sum_wcv += w[i]*g.Cv(Q)*g.M()[i] # where c_v is in J/mol.K
        
    dTdt = -sum_ew/sum_wcv
    dpdt = PC_R_u*Q.T[0]*sum(dwdt) + PC_R_u*dTdt*sum(w)
    dYdt = (array(dwdt)*g.M())/Q.rho
    
    if testFlag:
        print "printing time:"
        print t
        print "printing dTdt:"
        print "%10.9e" % dTdt
        print "printing dpdt:"
        print "%10.9e" % dpdt
        print "printing dYdt:"
        for dYidt in dYdt: 
            print "%10.9e" % dYidt
        sys.exit(1)
    dydt = array([dTdt] + [dpdt] + list(dYdt))

    return dydt

# ---------------------------------------------------------

def const_p_fixed_mass(y, t, Q, g, r, dt_chem, testFlag):
    # Y - mass fractions (kg/kg)
    # w - concentration (mol/m**3)
    # dwdt - molar production rate (mol/s)
    # M - molecular weights (kg/mol)
    
    nsp = g.get_number_of_species()

    # unpack array
    Q.T[0] = y[0] # temperature
    Q.p = y[1] # pressure
    Y = y[2:nsp+2] # mass fraction

    dYdt = zeros(nsp)

    dwdt = vectord(nsp*[0.0])
    w = vectord(nsp*[0.0])
    
    # update gas
    for isp in range(nsp):
        Q.massf[isp] = Y[isp]

    g.eval_thermo_state_pT(Q) # concs _are_ updated here
    
    # copy concs to w
    w = convert_massf2conc(Q.rho, Q.massf, g.M())

    # get rates-of-change
    dwdt = r.rate_of_change_py(Q)

    sum_hw = 0.0
    sum_wcp = 0.0

    for i in range(nsp):
        h_i = g.enthalpy(Q, i)*g.M()[i]
        sum_hw += h_i*dwdt[i]
        sum_wcp += w[i]*g.Cp(Q)*g.M()[i] # J/mol.K
        
    dTdt = -sum_hw/sum_wcp
    dpdt = 0.0
    dYdt = (array(dwdt)*g.M())/Q.rho

    if testFlag:
        print "printing time:"
        print t
        print "printing dTdt:"
        print "%10.9e" % dTdt
        print "printing dpdt:"
        print "%10.9e" % dpdt
        print "printing dYdt:"
        for dYidt in dYdt: 
            print "%10.9e" % dYidt
        sys.exit(1)
    dydt = array([dTdt] + [dpdt] + list(dYdt))

    return dydt

# ----------------------------------------------------------

def main():
    # handle input arguments
    userOptions, progArgs = getopt(sys.argv[1:], shortOptions, longOptions)

    if len(userOptions) == 0:
        printUsage()
        sys.exit(0)
    else:
        uoDict = dict(userOptions)

    if uoDict.has_key('--help') or uoDict.has_key('-h'):
        printUsage()
        sys.exit(0)

    typeFlag = uoDict.get('--type', None)
    chem = uoDict.get('--chem', None)
    gas = uoDict.get('--gas', None)
    mf = uoDict.get('--mf', None)
    X = uoDict.get('--X', None)
    T = uoDict.get('--T', None)
    p = uoDict.get('--p', None)
    rho = uoDict.get('--rho', None)
    tmax = uoDict.get('--tmax', None)
    dt = uoDict.get('--dt', None)

    # --------------------------------------------------
    
    ignFlag = False
    testFlag = False
    
    if typeFlag == 'cvfm': 
        test = const_v_fixed_mass
    elif typeFlag == 'cpfm': 
        test = const_p_fixed_mass
    else:
        print "Please define reactor type."
        print " --type={'cvfm'|'cpfm'}"
        sys.exit() 

    if chem == None:
        print "Please define mechanism, e.g."
        print " --chem=grimech30"
        sys.exit() 

    if mf == None and X == None:
        print "Please define mass or mole fractions, e.g."
        print " --mf=\"{'CH4', 1.0}\" or --X=..."
        sys.exit() 

    if tmax == None:
        # take one step
        dt = 1e-7
        tmax = 1e-7
        testFlag = True
    elif str(tmax) == 'tig':
        # stop at ignition
        ignFlag = True
        tmax = 2e-3
        if dt == None:
            dt = 1e-9
        else: dt = float(dt)
    else:
        # user specified
        tmax = float(tmax)
        if dt == None:
            dt = tmax/2000.
        else: dt = float(dt)
    
    # ------------------------------------------------------
    # define gas and state

    g = create_gas_model(gas)
    Q = Gas_data(g)

    if mf != None:
        mf = eval(mf) # convert to dictionary
        set_massf(Q, g, mf)
    if X != None:
        X = eval(X)
        set_molef(Q, g, X)

    #for i in range(len(g.M())):
    #    print(g.species_name(i), Q.massf[i])
    
    if p == None: # we use rho, p
        M = g.mixture_molecular_weight(Q)
        Q.rho = float(rho)*M
        Q.T[0] = float(T)
        g.eval_thermo_state_rhoT(Q)
    elif rho == None: # we use T, p
        Q.T[0] = float(T)
        Q.p = float(p)*p_atm
        g.eval_thermo_state_pT(Q)
    else:
        print("insufficient data")
        print("please supply T and {p, rho}")
        sys.exit(1)
    r = create_Reaction_update(chem, g)
    
    # ------------------------------------------------------
    # solve system and print

    year, month, day, hrs, mins, secs, a, a, a = localtime()
    
    print "# %s reactor using %s" % (typeFlag, chem)
    # status = os.system("echo '# bzr:' `bzr log|head -n2|tail -n1`")
    print "# %02i:%02i, %02i/%02i/%4i" % (hrs, mins, day, month, year)
    print "#"

    # call test with arguments
    print "# t, T, p",
    for i in range(len(g.M())):
        print(", '%s'" % g.species_name(i)),
    print
    
    # ------------------------------------------------------
    # using scipy ode-update
    
    t0 = arange(0.0, tmax+dt, dt)
    y0 = [Q.T[0], Q.p]
    for i in range(len(g.M())):
        y0.append(Q.massf[i]) # mass fractions

    args = (Q, g, r, dt, testFlag)
    y = odeint(test, y0, t0, args, rtol=1e-6, atol=1e-15, full_output=True)
    
    for i in range(len(y[0])):
        T = y[0][i][0]
        p = y[0][i][1]
        c = convert_massf2conc(Q.rho, y[0][i][2:], g.M())
        c_tot = sum(c)
        print "%.5e %.5e %.5e " % (t0[i], T, p),
        # convert to mol fractions
        for j in range(len(c)):
            print "%.5e " % (c[j]/c_tot),
        print

    print '#',y[1]['message']

# ----------------------------------------------------------

shortOptions = 'h'
longOptions = ['help', 'type=', 'chem=', 'gas=', 'mf=', 'X=', 'T=', 'rho=', 'p=', 'tmax=', 'dt=']

def printUsage():
    print "Performs a constant-volume fixed-mass reactor test."
    print "-or- a constant-pressure fixed-mass reactor test."
    print "Writes data to stdout."
    print
    print "Usage: ./reactor.py [options]"
    print " -h, --help         prints this help"
    print
    print "[options]"
    print " --type             'cvfm' or 'cpfm'"
    print " --chem             suffix for the chemistry data files"
    print " --gas              thermodynamic data file"
    print " --mf, X            mass, mole fractions e.g. \"{'CH4', 1.0}\""
    print " --T                temperature, K"
    print " --p --rho          pressure, atm or density, mol/m^3"
    print " --tmax             final time, s"
    print " --dt               timestep, s"
    print
    print "if tmax, dt are not set, only one step is taken."

    return

# ----------------------------------------------------------

if __name__ == '__main__':
    main()
