#! /usr/bin/env python

## \file    simple_reactor.py
## \author  Brendan O'Flaherty
## 
## \version 10 Jun 2008 -- File created
## 
## \brief   A reactor that uses split operations.
## ---------------------------------------------------------

import os
import sys
import ConfigParser
from getopt import getopt
from time import localtime

from numpy import array, arange, zeros

# from zero_finding import bisection_root
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *

# ---------------------------------------------------------

def get_absolute_enthalpy(Q, T, X, h0f):
    """
    parameters:
    Q   is a Gas_data object
    T   is the temperature we are evaluating h0 at
    X   is an array of species mole fractions and
    h0f are their formation enthalpies
    
    """
    
    h_0 = calculate_enthalpy_at_T(Q, g, T_atm)
    h_1 = calculate_enthalpy_at_T(Q, g, T)
    dh = (h_1 - h_0)
    M = mixture_molecular_weight(Q)
    
    h0f_mix = Q.c_tot*sum(X*h0f)

    h0 = dh + h0f_mix
    
    #print "h0 = dh + h0f_mix"
    #print "   = %.4e + %.4e" % (dh, h0f_mix)
    #print "   = %.4e" % (h0)
    #print
    
    return h0

def const_pT(dt, args):
    Q, g, r, dt_chem, ignFlag, p_args = args
    
    #perform chemical substeps
    dt_chem = r.update_state_py(Q, dt, dt_chem, g)

    # update state
    g.eval_thermo_state_pT(Q)
    args = [Q, g, r, dt_chem, ignFlag, p_args]
    return args

# ---------------------------------------------------------

def const_v_fixed_mass(dt, args):
    Q, g, r, dt_chem, ignFlag, p_args = args
    
    # update assuming const volume
    dt_chem = r.update_state_py(Q, dt, dt_chem, g)
    
    # update state
    g.eval_thermo_state_rhoe(Q)
    
    if ignFlag:
        i_CH4 = g.get_isp_from_species_name("CH4");
        molef = convert_massf2molef(Q.massf, g.M())
        if molef[i_CH4] < 1e-6:
            #ignition
            ignFlag = False
            Q.T[0] = Q.T[0] - 2000 + 298.15
            g.eval_thermo_state_pT(Q)
            
        #p_0, dp_max, p_old = p_args
        #p_old.insert(0, Q.p); p_old.pop()
        #dp = (-p_old[2] + 4*p_old[1] - 3*p_old[0])/(2*dt)
        #d2p = (p_old[2] - 2*p_old[1] + p_old[0])/(dt*dt)
        
        #if (dp > dp_max):
        #    dp_max = dp
        #elif (Q.p > p_0 and d2p < 0.0):
        #    # ignition
        #    return False
        #p_args = [p_0, dp_max, p_old]

    args = [Q, g, r, dt_chem, ignFlag, p_args]
    return args

# ---------------------------------------------------------

def get_reaction_rate(args):
    Q, g, r, dt_chem, ignFlag, p_args = args
    
    return args

# ---------------------------------------------------------

def const_p_fixed_mass(dt, args):
    Q, g, r, dt_chem, ignFlag, p_args = args

    p_0 = Q.p
    v_0 = (1.0/Q.rho)
    gamma = g.gamma(Q)
    
    # update assuming const volume
    dt_chem = r.update_state_py(Q, dt, dt_chem, g)
    g.eval_thermo_state_rhoe(Q)

    # expand isentropically
    dp = Q.p - p_0
    v_1 = v_0*(Q.p/p_0)**(1.0/gamma)
    dv = v_1 - v_0
    de = 0.5*(p_0 + Q.p)*dv
    
    Q.e[0] -= de
    Q.rho *= v_0/v_1

    # update state
    g.eval_thermo_state_rhoe(Q)

    #print("# de=%12.11e, dv=%12.11e, dp=%12.11e" % (dv, de, dp))
    #print("# Q.p=%12.11e p_0=%12.11e" % (Q.p, p_0))

    # -----------------------------------------------------

    if ignFlag:
        p_0, dp_max, p_old = p_args

        p_old.insert(0, Q.p); p_old.pop()
        dp = (-p_old[2] + 4*p_old[1] - 3*p_old[0])/(2*dt)
        d2p = (p_old[2] - 2*p_old[1] + p_old[0])/(dt*dt)
        
        if (dp > dp_max):
            dp_max = dp
        elif (Q.p > p_0 and d2p < 0.0):
            # ignition has begun
            return False
        p_args = [p_0, dp_max, p_old]

    args = [Q, g, r, dt_chem, ignFlag, p_args]
    return args

# ----------------------------------------------------------

def run_reactor(Q, g, r, dt, tmax, reactor, testFlag, ignFlag, writeFlag=True):
    # set initial vectors

    t0 = arange(0.0, tmax, dt)
    dt_chem = -1.0
    T0 = Q.T[0]
    p0 = Q.p

    # ------------------------------------------------------
    # work done here
    print "# t, T_0, p_0, T, p, rho",
    for i in range(len(g.M())):
        print(", '%s'" % g.species_name(i)),
    print

    #count = 0
    if writeFlag: write_state(t0[0], Q, g, T0, p0)
    if ignFlag:
        p_old = [Q.p, Q.p, Q.p]
        p_args = [Q.p, 0.0, p_old]
        args = [Q, g, r, dt_chem, ignFlag, p_args]
        
        for t in t0[1:]:
            args = reactor(dt, args)
            if type(args) == bool:
                write_state(t, Q, g, T0, p0)
                write_gas_state(t, Q, g)
                break
            if writeFlag: write_state(t, Q, g, T0, p0)
        if t == t0[-1]: 
            write_gas_state(t, Q, g)
            print  "# ignition not reached by t = ", t
    else: 
        args = [Q, g, r, dt_chem, ignFlag, []]
        for t in t0[1:]:
            args = reactor(dt, args)
            if writeFlag: 
                Q, g, r, dt_chem, ignFlag, p_args = args
                write_state(t, Q, g, T0, p0)
    return

# ----------------------------------------------------------

def write_state(t, Q, g, T0, p0):
    print("%.5e %.5e %.5e %.5e %.5e %.12e" % (t, T0, p0, Q.T[0], Q.p, Q.rho)),
    molef = convert_massf2molef(Q.massf, g.M())
    for i in range(len(g.M())):
        print(" %.12e" % molef[i]),
    print
    return

def write_gas_state(t, Q, g):
    molef = convert_massf2molef(Q.massf, g.M())
    print("--X=\"{"),
    for i in range(len(g.M())):
        print("'%s':%.12e," % (g.species_name(i), molef[i])),
    print("}\"")
    return

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
    ignFlag = uoDict.get('--ign', False)
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

    testFlag = False

    if typeFlag == None:
        print "Please define reactor type."
        print " --type={'cvfm'|'cpfm'|'ctp'}"
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
    
    # -----------------------------------------------------
    # define gas and state

    g = create_gas_model(gas)
    Q = Gas_data(g)
    # not needed any more # g.initialise_gas_data(Q)

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
    
    # -----------------------------------------------------
    # select reactor

    if typeFlag == 'ctp': 
        print "# const_pT reactor using %s" % (chem)
        reactor = const_pT
    elif typeFlag == 'cvfm': 
        print "# const_v_fixed_mass reactor using %s" % (chem)
        reactor = const_v_fixed_mass
    elif typeFlag == 'cpfm': 
        print "# const_p_fixed_mass reactor using %s" % (chem)
        if ignFlag:
            print "ignition is calculated using POI of pressure."
            print "it makes no sense to run a cpfm reactor to ignition."
            sys.exit(1)
        reactor = const_p_fixed_mass
    else:
        print "Reactor type", typeFlag, "unknown."
        print "Exiting to system..."
        sys.exit(1)

    # -----------------------------------------------------
    # print preamble
    
    year, month, day, hrs, mins, secs, a, a, a = localtime()
    
    #status = os.system("echo '# git' `git log | head -n1`")
    print "# %02i:%02i, %02i/%02i/%4i" % (hrs, mins, day, month, year)
    print "#"

    # -----------------------------------------------------
    # get adiabatic flame temperature

#    # CH4, H2O, CO2, O2, N2
#    h0f = array([-74831., -241845., -393546., 0., 0.]) # J/mol
#    X_r = array([1.0, 0.0, 0.0, 2.0, 7.52])/10.52
#    X_p = array([0.0, 2.0, 1.0, 0.0, 7.52])/10.52
#    
#    # create burnt gas
#    Q_p = Gas_data()
#    moles_p = [['CO2', 1.0], ['H2O', 2.0], ['N2', 7.52]]
#    Q_p, mix_p = define_mixture(Q_p, species_list, moles_p)
#    Q_p.T = Q.T[0]
#    Q_p.p = Q.p
#    g.eval_thermo_state_pT(Q_p)
#
#    # calculate absolute enthalpy of reactants
#    #h0_r = get_absolute_enthalpy(Q, T, X_r, h0f)
#    h0_r = Q.rho*calculate_absolute_enthalpy(Q)
#
#    args = (Q_p, X_p, h0f, h0_r)
#    
#    # iterate to find absolute enthalpy of products
#    print "T0 = %.4e" % (Q.T[0])
#    Tf = bisection_root(get_dh0, Q.T[0], Q.T[0]+2e3, tol=1e-3, args=args)
#    print "Tf = %.4e" % (Tf)
    
    # -----------------------------------------------------
    # call test with arguments
    run_reactor(Q, g, r, dt, tmax, reactor, testFlag, ignFlag)
    print "# Done."

# ----------------------------------------------------------

def get_dh0(T_guess, args):
    Q_p, X_p, h0f, h0_r = args
    #h0_p = get_absolute_enthalpy(Q_p, T_guess, X_p, h0f)
    h0_p = Q_p.rho*calculate_absolute_enthalpy(Q_p, T_guess)

    return h0_p - h0_r

# ----------------------------------------------------------

shortOptions = 'h'
longOptions = ['help', 'type=', 'chem=', 'ign=', 'gas=', 'mf=', 'X=', 'T=', 'rho=', 'p=', 'tmax=', 'dt=']
    
def printUsage():
    print "Performs a constant-volume fixed-mass reactor test."
    print "-or- a constant-pressure fixed-mass reactor test."
    print "Writes data to stdout."
    print
    print "Usage: ./simple_reactor.py [options]"
    print " -h, --help        prints this help"
    print
    print "[options]"
    print " --type             'cvfm', 'cpfm' or 'ctp'"
    print " --chem             suffix for the chemistry data files"
    print " --gas              thermodynamic data file"
    print " --ign              run the reaction until ignition or tmax, whichever occurs first"
    print " --mf, X            mass, mole fractions e.g. \"{'CH4', 1.0}\""
    print " --T                temperature, K"
    print " --p --rho          pressure, atm or density, mol/m^3"
    print " --tmax             final time, s"
    print " --dt               timestep, s"
    print
    print "if tmax, dt are not set, only one step is taken."
    print "if tmax='tig', then program breaks at maximum dpdt."

    return

# ----------------------------------------------------------

if __name__ == '__main__':
    main()
