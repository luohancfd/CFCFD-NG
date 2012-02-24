#! /usr/bin/env python 
"""
estcj.py: Equilibrium Shock Tube Conditions, Junior

Since 1968, we have been using the ESTC code by Malcolm McIntosh
to compute the conditions in the end of the reflected shock tubes
T1--T5 and HEG.  There are a number of problems in using the ESTC
code, including uncertainty in updating the chemistry coefficients.

This program, ESTCJ, moves away from the old chemistry model
by making use of the CEA code from the NASA Glenn Research Center.

.. Author: PA Jacobs
   Institute of Aerodynamics and Flow Technology
   The German Aerospace Center, Goettingen.

.. Versions:
   24-Dec-02: First code.
   2010: ported to run with Rowan's cea2_gas module.
   2011: Added isentropic expansions so that we now have
       a full replacement for stn.f
   01-June-2011 LukeD: Separated the code which writes an output
       file into its own function to allow for better integration with nenzfr.py
   30-June-2011 LukeD: Decreased the starting guess for secant 
       when solving for the exit flow
   22-July-2011 LukeD: Added stnp option which allows us to expand
       to a nominated pitot-to-supply pressure ratio. The calculated pitot 
       pressure and pitot-to-supply pressure ratio are included in the values 
       printed out for the nozzle exit
   24-Feb-2012 PJ: update to use the new cea2_gas.py arrangement.
"""

import sys, os, math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from cfpylib.nm.zero_solvers import secant
# We base our calculation of gas properties upon calls to the NASA Glenn CEA code.
from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.ideal_gas_flow import *

# ----------------------------------------------------------------------------

VERSION_STRING = "24-Feb-2012"
DEBUG_ESTCJ  = 0  # if 1: some detailed data is output to help debugging
PRINT_STATUS = 1  # if 1: the start of each stage of the computation is noted.

# ----------------------------------------------------------------------------
# Utility functions.

def make_gas_from_name(gasName):
    """
    Manufacture a cea2_gas object from a small library of options.

    :param gasName: one of the names for the special cases set out below
    """
    if gasName == 'air':
        return Gas({'Air':1.0,})
    elif gasName == 'air5species':
        return Gas(reactants={'N2':0.79, 'O2':0.21, 'N':0.0, 'O':0.0, 'NO':0.0}, 
                   inputUnits='moles', onlyList=['N2','O2','N','O','NO'])
    elif gasName == 'n2':
        return Gas(reactants={'N2':1.0, 'N':0.0}, onlyList=['N2', 'N'])
    elif gasName == 'co2':
        return Gas(reactants={'CO2':1.0})
    elif gasName == 'h2ne':
        return Gas(reactants={'H2':0.15, 'Ne':0.85}, inputUnits='moles')
    else:
        raise Exception, 'make_gas_from_name(): unknown gasName: %s' % gasName


def shock_ideal(s1, Vs, s2):
    """
    Computes post-shock conditions in the shock frame, assuming ideal gas.

    :param s1: pre-shock gas state
    :param Vs: speed of gas coming into shock
    :param s2: post-shock gas state
    :returns: both the shock-reference speed and the lab speed.
    """
    #
    M1 = Vs / s1.a
    V1 = Vs
    gam = s1.gam
    R = s1.R
    C_v = s1.C_v
    #
    s2.rho = s1.rho * (gam + 1.0) * M1 * M1 \
             / (2.0 + (gam - 1.0) * M1 * M1)
    s2.p = s1.p * (2.0 * gam * M1 * M1 - (gam - 1.0)) / (gam + 1.0)
    s2.T = s2.p / (R * s2.rho)
    s2.e = s2.T * C_v
    #
    V2 = s1.rho / s2.rho * V1
    Vg = V1 - V2
    s2.a = s1.a * math.sqrt(s2.T / s1.T)
    #
    s2.R = s1.R
    s2.gam = s1.gam
    s2.C_v = s1.C_v
    #
    return (V2,Vg)

def my_limiter(delta, orig, frac=0.5):
    """
    Limit the magnitude of delta to no more than a fraction of the original.
    """
    if delta >= 0.0:
        sign = 1
    else:
        sign = -1
    abs_delta = min(abs(delta), frac*abs(orig))
    return sign * abs_delta

def shock_real(s1, Vs, s2):
    """
    Computes post-shock conditions, using high-temperature gas properties
    and a shock-stationary frame.

    :param s1: pre-shock gas state
    :param Vs: speed of gas coming into shock
    :param s2: post-shock gas state
    :returns: both the shock-reference speed and the lab speed.
    """
    #
    (V2,Vg) = shock_ideal(s1, Vs, s2)
    if DEBUG_ESTCJ:
        print 'shock_real(): post-shock condition assuming ideal gas'
        s2.write_state(sys.stdout)
        print '    V2: %g m/s, Vg: %g m/s' % (V2,Vg)
    #
    # We assume that p1 and T1 are correct
    # and that s2 contains a fair initial guess.
    V1 = Vs
    s1.EOS(problemType='pT');
    if DEBUG_ESTCJ:
        print 'shock_real(): pre-shock condition assuming real gas and original pT'
        s1.write_state(sys.stdout)
    s2.EOS(problemType='pT');
    if DEBUG_ESTCJ:
        print 'shock_real(): post-shock condition assuming real gas and ideal pT'
        s2.write_state(sys.stdout)
    #
    p1    = s1.p
    e1    = s1.e
    rho1  = s1.rho
    temp1 = p1 + rho1 * V1 * V1;             # momentum
    H1    = e1 + p1 / rho1 + 0.5 * V1 * V1;  # total enthalpy
    #
    rho_delta = 1.0
    T_delta = 1.0
    rho_tol = 1.0e-3; # tolerance in kg/m^3
    T_tol = 0.25;  # tolerance in degrees K
    #
    s0 = s1.clone(); # temporary workspace
    #
    # Update the estimates using the Newton-Raphson method.
    #
    for count in range(20):
        rho_save = s2.rho
        T_save = s2.T
        p_save = s2.p
        #
        p2 = p_save
        T2 = T_save
        s2.EOS(problemType='pT')
        e2 = s2.e
        rho2 = rho_save
        r2r1 = rho2 / rho1
        f1 = temp1 - p2 - rho1 * rho1 * V1 * V1 / rho2
        f2 = H1 - e2 - p2 / rho2 - 0.5 * V1 * V1 / (r2r1 * r2r1)
        f1_save = f1
        f2_save = f2
        #
        # Use finite differences to compute the Jacobian.
        d_rho = rho_save * 0.01
        d_T = T_save * 0.01
        #
        rho2 = rho_save + d_rho
        T2 = T_save
        s0.rho = rho2
        s0.T = T2
        s0.EOS(problemType='rhoT')
        p2 = s0.p
        e2 = s0.e
        f1 = temp1 - p2 - rho1 * rho1 * V1 * V1 / rho2
        f2 = H1 - e2 - p2 / rho2 - 0.5 * V1 * V1 / (r2r1 * r2r1)
        A = (f1 - f1_save) / d_rho
        C = (f2 - f2_save) / d_rho
        #
        rho2 = rho_save
        T2 = T_save + d_T
        s0.rho = rho2
        s0.T = T2
        s0.EOS(problemType='rhoT')
        e2 = s0.e
        p2 = s0.p
        f1 = temp1 - p2 - rho1 * rho1 * V1 * V1 / rho2
        f2 = H1 - e2 - p2 / rho2 - 0.5 * V1 * V1 / (r2r1 * r2r1)
        B = (f1 - f1_save) / d_T
        D = (f2 - f2_save) / d_T
        #
        # Invert Jacobian and multiply.
        det = A * D - B * C
        rho_delta = (D * f1_save - B * f2_save) / det
        T_delta = (-C * f1_save + A * f2_save) / det
        #
        rho_delta = my_limiter(rho_delta, rho_save)
        T_delta = my_limiter(T_delta, T_save)
        rho_new = rho_save - rho_delta
        T_new   = T_save - T_delta
        if DEBUG_ESTCJ:
            print('shock_real(): rho_save=%e, T_save=%e' % (rho_save, T_save))
            print('shock_real(): rho_delta=%e, T_delta=%e' % (rho_delta, T_delta))
            print('shock_real(): rho_new=%e, T_new=%e' % (rho_new, T_new))
        #
        s2.rho = rho_new
        s2.T   = T_new
        s2.EOS(problemType='rhoT')
        #
        # Check convergence.
        if abs(rho_delta) < rho_tol and abs(T_delta) < T_tol: break
    #
    if DEBUG_ESTCJ:
        print ('shock_real(): count = %d, drho=%e, dT=%e' %
               (count, rho_delta, T_delta) )
    #
    # Back-out velocities via continuity.
    V2 = V1 * s1.rho / s2.rho
    Vg = V1 - V2
    return (V2, Vg)


def reflected_shock(s2, Vg, s5):
    """
    Computes state5 which has brought the gas to rest at the end of the shock tube.

    :param state2: the post-incident-shock gas state
    :param Vg: the lab-frame velocity of the gas in state 2
    :param s5: the stagnation state that will be filled in
        (as a side effect of this function)
    :returns: Vr, the reflected shock speed in the lab frame.
    """
    #
    # As an initial guess, 
    # assume that we have a very strong shock in an ideal gas.
    density_ratio = (s2.gam + 1.0)/(s2.gam - 1.0)
    Vr_a = Vg / density_ratio;
    V5, Vjunk = shock_real(s2, Vr_a+Vg, s5)
    # The objective function is the difference in speeds,
    # units are m/s.  A value of zero for this function means
    # that, as the shock propagates upstream with speed ur,
    # the processed test gas is left in the end of the tube
    # with a velocity of zero in the laboratory frame.
    f_a = V5 - Vr_a
    if DEBUG_ESTCJ:
        print 'Reflected shock: Vr_a: %g, V5: %g' % (Vr_a, V5)
    #
    # Now, we need to update this guess...use a secant update.
    #
    Vr_b = 1.1 * Vr_a
    V5, Vjunk = shock_real(s2, Vr_b+Vg, s5)
    f_b = V5 - Vr_b
    if DEBUG_ESTCJ:
        print 'Reflected shock: Vr_b: %g, V5: %g' % (Vr_b, V5)
    if abs(f_a) < abs(f_b):
        f_a, f_b = f_b, f_a
        Vr_a, Vr_b = Vr_b, Vr_a
    count = 0
    while abs(f_b) > 0.5 and count < 20:
        slope = (f_b - f_a) / (Vr_b - Vr_a)
        Vr_c = Vr_b - f_b / slope
        V5, Vjunk = shock_real(s2, Vr_c+Vg, s5)
        f_c = V5 - Vr_c
        if abs(f_c) < abs(f_b):
            Vr_b = Vr_c; f_b = f_c
        else:
            Vr_a = Vr_c; f_a = f_c
        count = count + 1
    #
    # At this point, ur_b should be out best guess.
    # Update the gas state data and return the best-guess value.
    #
    if count >= 20:
        print 'Reflected shock iteration did not converge.'
    V5, Vjunk = shock_real(s2, Vr_b+Vg, s5)
    return Vr_b


def expand_from_stagnation(p_over_p0, p0, T0, gasName):
    """
    Given a stagnation condition, expand to a new pressure.
    Returns new gas state and the corresponding velocity of the expanded stream.
    """
    state0 = make_gas_from_name(gasName)
    state0.set_from_pAndT(p0, T0)
    new_state = make_gas_from_name(gasName)
    new_state.set_from_pAndT(p0, T0)
    new_state.p = p0 * p_over_p0;
    new_state.EOS(problemType='ps')
    # Matt McGilvray had a note about CEA giving bad entropy values
    # so we'll assert things are OK before proceeding.
    assert abs(new_state.s - state0.s)/abs(state0.s) < 0.001
    h = new_state.e + new_state.p/new_state.rho  # static enthalpy
    H = state0.e + state0.p/state0.rho  # stagnation enthalpy
    V = math.sqrt(2.0*(H-h))
    return new_state, V


def total_condition(p1, T1, V1, gasName):
    """
    Given a free-stream condition, compute the corresponding stagnant condition
    at which the gas is brought to rest isentropically.
    """
    state1 = make_gas_from_name(gasName)
    state1.set_pT(p1, T1)
    H1 = state1.p/state1.rho + state1.e + 0.5*V1*V1
    def error_in_total_enthalpy(x, p1=state1.p, T1=state1.T, s1=state1.s, H1=H1):
        """
        The enthalpy at the stagnation condition should match
        the total enthalpy of the stream.
        """
        new_state = make_gas_from_name(gasName)
        new_state.set_pT(p1, T1)
        new_state.p = x * p1
        new_state.EOS(problemType='ps')
        h = new_state.p/new_state.rho + new_state.e
        return (H1 - h)/abs(H1)
    x_total = secant(error_in_total_enthalpy, 1.0, 1.01, tol=1.0e-4)
    if x_total == 'FAIL':
        print "Failed to find total conditions iteratively."
        x_total = 1.0
    new_state = make_gas_from_name(gasName)
    new_state.set_pT(p1, T1)
    new_state.p = x_total * p1
    new_state.EOS(problemType='ps')
    return new_state


def pitot_condition(p1, T1, V1, gasName):
    """
    Given a free-stream condition, compute the corresponding Pitot condition
    at which the gas is brought to rest, possibly through a shock.
    """
    state1 = make_gas_from_name(gasName)
    state1.set_pT(p1, T1)
    if V1 > state1.a:
        # Supersonic free-stream; process through a shock first.
        state2 = make_gas_from_name(gasName)
        (V2,Vg) = shock_real(state1, V1, state2)
        return total_condition(state2.p, state2.T, V2, gasName)
    else:
        # Subsonic free-stream
        return total_condition(p1, T1, V1, gasName)
    
#--------------------------------------------------------------------

def reflected_shock_tube_calculation(gasName, p1, T1, Vs, pe, pp_on_pe, area_ratio, task):
    """
    Runs the reflected-shock-tube calculation from initial fill conditions
    observed shock speed and equilibrium pressure.

    :param task: one of 'ishock', 'st', 'stn', 'stnp'
    """
    if PRINT_STATUS: print 'Write pre-shock condition.'
    state1 = make_gas_from_name(gasName)
    state1.set_pT(p1, T1)
    H1 = state1.u + state1.p/state1.rho
    result = {'state1':state1, 'H1':H1}
    #
    if PRINT_STATUS: print 'Start incident-shock calculation.'
    state2 = make_gas_from_name(gasName)
    (V2,Vg) = shock_real( state1, Vs, state2 )
    result['state2'] = state2
    result['V2'] = V2
    result['Vg'] = Vg
    #
    if task == 'ishock':
        # We want post-incident-shock conditions only.
        return result
    #
    if PRINT_STATUS: print 'Start reflected-shock calculation.'
    state5 = make_gas_from_name(gasName)
    Vr = reflected_shock(state2, Vg, state5)
    result['state5'] = state5
    result['Vr'] = Vr
    #
    if PRINT_STATUS: print 'Start calculation of isentropic relaxation.'
    state5s = make_gas_from_name(gasName)
    state5s.set_pT(state5.p, state5.T);  # entropy is set
    state5s.p = pe;                    # then pressure is relaxed
    state5s.EOS(problemType='ps');     # via an isentropic process
    result['state5s'] = state5s
    H5s = state5s.e + state5s.p/state5s.rho # stagnation enthalpy
    result['H5s'] = H5s
    #
    if task in ['stn','stnp']:
        if PRINT_STATUS: print 'Start isentropic relaxation to throat (Mach 1)'
        def error_at_throat(x, s5s=state5s, gasName=gasName):
            "Returns Mach number error as pressure is changed."
            state, V = expand_from_stagnation(x, s5s.p, s5s.T, gasName)
            return (V/state.a) - 1.0
        x6 = secant(error_at_throat, 0.95, 0.90, tol=1.0e-4)
        if x6 == 'FAIL':
            print "Failed to find throat conditions iteratively."
            x6 = 1.0
        state6, V6 = expand_from_stagnation(x6, state5s.p, state5s.T, gasName)
        mflux6 = state6.rho * V6  # mass flux per unit area, at throat
        result['state6'] = state6
        result['V6'] = V6
        result['mflux6'] = mflux6
        #
        if task == 'stn':
            if PRINT_STATUS: print 'Start isentropic relaxation to nozzle exit.'
            # The mass flux going through the nozzle exit has to be the same
            # as that going through the nozzle throat.
            def error_at_exit(x, s5s=state5s, s6=state6, mflux_throat=mflux6,
                           area_ratio=area_ratio, gasName=gasName):
                "Returns mass_flux error as pressure is changed."
                state, V = expand_from_stagnation(x, s5s.p, s5s.T, gasName)
                mflux = state.rho * V * area_ratio
                if DEBUG_ESTCJ: print "x=", x, "p=", state.p, "T=", state.T, "V=", V, \
                        "mflux=", mflux, "mflux_throat=", mflux_throat
                return (mflux-mflux_throat)/mflux_throat
            # It appears that we need a pretty good starting guess for the pressure ratio.
            # Maybe a low value is OK.
            x7 = secant(error_at_exit, 0.001*x6, 0.00005*x6, tol=1.0e-4)
            if x7 == 'FAIL':
                print "Failed to find exit conditions iteratively."
                x7 = x6
            state7, V7 = expand_from_stagnation(x7, state5s.p, state5s.T, gasName)
            mflux7 = state7.rho * V7 * area_ratio
            result['area_ratio'] = area_ratio
            state7_pitot = pitot_condition(state7.p, state7.T, V7, gasName)
            result['state7'] = state7
            result['V7'] = V7
            result['mflux7'] = mflux7
            result['pitot7'] = state7_pitot.p
        elif task == 'stnp':
            if PRINT_STATUS: print 'Start isentropic relaxation to nozzle exit pitot pressure.'
            # The exit pitot pressure has to be the same as that measured
            def error_at_exit(x, s5s=state5s, s6=state6, pp_pe=pp_on_pe, gasName=gasName):
                "Returns pitot pressure error as static pressure is changed."
                state1, V = expand_from_stagnation(x, s5s.p, s5s.T, gasName)
                state2 = pitot_condition(state1.p, state1.T, V, gasName)
                if DEBUG_ESTCJ: print "x=", x, "pitot_to_supply=", state2.p/s5s.p, \
                    "relative error=", (state2.p/s5s.p - pp_pe)/pp_pe
                return (state2.p/s5s.p - pp_pe)/pp_pe
            # We need a low starting guess for the pressure ratio.
            #x7 = secant(error_at_exit, 0.001*x6, 0.00005*x6, tol=1.0e-4)
            # Changed the tolerance on 25/07/2011 in order to get the M8 nozzle to work (shot 10803)
            x7 = secant(error_at_exit, 0.001*x6, 0.00005*x6, tol=1.5e-4)
            if x7 == 'FAIL':
                print "Failed to find exit conditions iteratively."
                x7 = x6
            state7, V7 = expand_from_stagnation(x7, state5s.p, state5s.T, gasName)
            result['area_ratio'] = mflux6/(state7.rho * V7)
            state7_pitot = pitot_condition(state7.p, state7.T, V7, gasName)
            #mflux7 = mflux6
            result['state7'] = state7
            result['V7'] = V7
            result['mflux7'] = mflux6
            result['pitot7'] = state7_pitot.p
            if DEBUG_ESTCJ: print "area_ratio=", area_ratio, "pitot7=", state7_pitot.p
    #
    if PRINT_STATUS: print 'Done with reflected shock tube calculation.'
    return result

#--------------------------------------------------------------------

def main():
    """
    The application gets information from the command options,
    does some calculation (depending on the specified task)
    and writes the results to the console or a file.
    """
    import optparse
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--task', dest='task', default='st',
                  choices=['st', 'stn', 'stnp', 'ishock', 'total', 'pitot'],
                  help=("particular calculation to make: "
                        "st = reflected shock tube; "
                        "stn = reflected shock tube with nozzle; "
                        "stnp = reflected shock tube with nozzle expanded to pitot; "
                        "ishock = incident shock only; "
                        "total = free-stream to total conditions; "
                        "pitot = free-stream to Pitot condition"))
    op.add_option('--gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'co2', 'h2ne'],
                  help=("name of gas model: "
                        "air; " "air5species; " "n2; " "co2; " "h2ne"))
    op.add_option('--p1', dest='p1', type='float', default=None,
                  help=("shock tube fill pressure or static pressure, in Pa"))
    op.add_option('--T1', dest='T1', type='float', default=None,
                  help=("shock tube fill temperature, in degrees K"))
    op.add_option('--V1', dest='V1', type='float', default=None,
                  help=("initial speed of gas in lab frame [default: %default], in m/s"))
    op.add_option('--Vs', dest='Vs', type='float', default=None,
                  help=("incident shock speed, in m/s"))
    op.add_option('--pe', dest='pe', type='float', default=None,
                  help=("equilibrium pressure (after shock reflection), in Pa"))
    op.add_option('--pp_on_pe', dest='pp_on_pe', type='float', default=None,
                  help=("nozzle supply to exit pitot pressure ratio"))
    op.add_option('--ar', dest='area_ratio', type='float', default=None,
                  help=("exit-to-throat area ratio of the nozzle"))
    op.add_option('--ofn', dest='outFileName', default=None,
                  help="name of file in which to accumulate output."
                      " file name will be: outFileName-estcj.dat"
                      " (Note that output defaults to stdout.)")
    opt, args = op.parse_args()
    #
    task = opt.task
    gasName = opt.gasName
    p1 = opt.p1
    T1 = opt.T1
    V1 = opt.V1
    Vs = opt.Vs
    pe = opt.pe
    pp_on_pe = opt.pp_on_pe
    area_ratio = opt.area_ratio
    outFileName = opt.outFileName
    if DEBUG_ESTCJ: print 'estcj:', gasName, p1, T1, V1, Vs, pe, area_ratio, outFileName
    #
    bad_input = False
    if p1 is None:
        print "Need to supply a float value for p1."
        bad_input = True
    if T1 is None:
        print "Need to supply a float value for T1."
        bad_input = True
    if Vs is None and task in ['stn', 'stnp', 'st', 'ishock']:
        print "Need to supply a float value for Vs."
        bad_input = True
    if pe is None and task in ['stn', 'stnp', 'st']:
        print "Need to supply a float value for pe."
        bad_input = True
    if pp_on_pe is None and task in ['stnp']:
        print "Need to supply a float value for pp_on_pe."
        bad_input = True
    if area_ratio is None and task in ['stn']:
        print "Need to supply a float value for ar=area_ratio."
        bad_input = True
    if bad_input:
        return -2
    #
    if outFileName is None:
        fout = sys.stdout
    else:
        fout = open(outFileName+'-estcj.dat','w')
    fout.write('estcj: Equilibrium Shock Tube Conditions\n')
    fout.write('Version: %s\n' % VERSION_STRING)
    #
    if task in ['st', 'stn', 'stnp', 'ishock']:
        fout.write('Input parameters:\n')
        fout.write('    Gas is %s, p1: %g Pa, T1: %g K, Vs: %g m/s\n' % (gasName,p1,T1,Vs) )
        result = reflected_shock_tube_calculation(gasName, p1, T1, Vs,
                                                  pe, pp_on_pe, area_ratio,
                                                  task=task)
        fout.write('State 1: pre-shock condition\n')
        result['state1'].write_state(fout)
        fout.write('State 2: post-shock condition.\n')
        result['state2'].write_state(fout)
        fout.write('    V2: %g m/s, Vg: %g m/s\n' % (result['V2'],result['Vg']) )
        if task in ['st', 'stn', 'stnp']:
            fout.write('State 5: reflected-shock condition.\n')
            result['state5'].write_state(fout)
            fout.write('    Vr: %g m/s\n' % (result['Vr'],) )
            fout.write('State 5s: equilibrium condition (relaxation to pe)\n')
            result['state5s'].write_state(fout)
            fout.write('Enthalpy difference (H5s - H1): %g J/kg\n' % 
                       ((result['H5s'] - result['H1']),) )
            if task in ['stn','stnp']:
                # shock tube plus nozzle, expand gas isentropically, stopping at area_ratio
                fout.write('State 6: Nozzle-throat condition (relaxation to M=1)\n')
                result['state6'].write_state(fout)
                fout.write('    V6: %g m/s, M6: %g, mflux6: %g kg/s/m**2\n' % 
                           (result['V6'], result['V6']/result['state6'].a, result['mflux6'],) )
                fout.write('State 7: Nozzle-exit condition (relaxation to correct mass flux)\n')
                result['state7'].write_state(fout)
                fout.write('    V7: %g m/s, M7: %g, mflux7: %g kg/s/m**2, area_ratio: %g, pitot: %g Pa\n' % 
                           (result['V7'], result['V7']/result['state7'].a, result['mflux7'],
                            result['area_ratio'], result['pitot7'],) )
                fout.write('    pitot7_on_p5s: %g\n' % (result['pitot7']/result['state5s'].p,) )
    elif task in ['total', 'TOTAL', 'Total']:
        fout.write('Input parameters:\n')
        fout.write('    Gas is %s, p1: %g Pa, T1: %g K, V1: %g m/s\n' % (gasName,p1,T1,V1) )
        state0 = total_condition(p1, T1, V1, gasName)
        fout.write('Total condition:\n')
        state0.write_state(fout)
    elif task in ['pitot', 'PITOT', 'Pitot']:
        fout.write('Input parameters:\n')
        fout.write('    Gas is %s, p1: %g Pa, T1: %g K, V1: %g m/s\n' % (gasName,p1,T1,V1) )
        state0 = pitot_condition(p1, T1, V1, gasName)
        fout.write('Pitot condition:\n')
        state0.write_state(fout)
    #
    if outFileName is None:
        pass
    else:
        fout.close()
    return 0

# -------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "Equilibrium Shock Tube Conditions"
        print "   Version:", VERSION_STRING
        print "   To see some useful hints, invoke this program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
