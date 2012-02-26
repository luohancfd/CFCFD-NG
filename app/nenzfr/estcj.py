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
from cfpylib.gasdyn.cea2_gas_flow import shock_real, reflected_shock
from cfpylib.gasdyn.cea2_gas_flow import expand_from_stagnation, total_condition, pitot_condition
from cfpylib.gasdyn.ideal_gas_flow import *

# ----------------------------------------------------------------------------

VERSION_STRING = "26-Feb-2012"
DEBUG_ESTCJ  = 0  # if 1: some detailed data is output to help debugging
PRINT_STATUS = 1  # if 1: the start of each stage of the computation is noted.

# ----------------------------------------------------------------------------

def make_gas_from_name(gasName, outputUnits='massf'):
    """
    Manufacture a cea2_gas object from a small library of options.

    :param gasName: one of the names for the special cases set out below
    """
    if gasName == 'air':
        return Gas({'Air':1.0,}, outputUnits=outputUnits)
    elif gasName == 'air5species':
        return Gas(reactants={'N2':0.79, 'O2':0.21, 'N':0.0, 'O':0.0, 'NO':0.0}, 
                   inputUnits='moles', onlyList=['N2','O2','N','O','NO'],
                   outputUnits=outputUnits)
    elif gasName == 'n2':
        return Gas(reactants={'N2':1.0, 'N':0.0}, onlyList=['N2', 'N'],
                   outputUnits=outputUnits)
    elif gasName == 'co2':
        return Gas(reactants={'CO2':1.0}, outputUnits=outputUnits)
    elif gasName == 'h2ne':
        return Gas(reactants={'H2':0.15, 'Ne':0.85}, inputUnits='moles',
                   outputUnits=outputUnits)
    else:
        raise Exception, 'make_gas_from_name(): unknown gasName: %s' % gasName

#--------------------------------------------------------------------

def reflected_shock_tube_calculation(gasName, p1, T1, Vs, pe, pp_on_pe, area_ratio, task):
    """
    Runs the reflected-shock-tube calculation from initial fill conditions
    observed shock speed and equilibrium pressure.

    This function may be imported into other applications (such as nenzfr).

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
    (V2,Vg) = shock_real(state1, Vs, state2)
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
            state, V = expand_from_stagnation(x, s5s)
            return (V/state.a) - 1.0
        x6 = secant(error_at_throat, 0.95, 0.90, tol=1.0e-4)
        if x6 == 'FAIL':
            print "Failed to find throat conditions iteratively."
            x6 = 1.0
        state6, V6 = expand_from_stagnation(x6, state5s)
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
                           area_ratio=area_ratio):
                "Returns mass_flux error as pressure is changed."
                state, V = expand_from_stagnation(x, s5s)
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
            state7, V7 = expand_from_stagnation(x7, state5s)
            mflux7 = state7.rho * V7 * area_ratio
            result['area_ratio'] = area_ratio
            state7_pitot = pitot_condition(state7, V7)
            result['state7'] = state7
            result['V7'] = V7
            result['mflux7'] = mflux7
            result['pitot7'] = state7_pitot.p
        elif task == 'stnp':
            if PRINT_STATUS: print 'Start isentropic relaxation to nozzle exit pitot pressure.'
            # The exit pitot pressure has to be the same as that measured
            def error_at_exit(x, s5s=state5s, s6=state6, pp_pe=pp_on_pe):
                "Returns pitot pressure error as static pressure is changed."
                state1, V = expand_from_stagnation(x, s5s)
                state2 = pitot_condition(state1, V)
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
            state7, V7 = expand_from_stagnation(x7, state5s)
            result['area_ratio'] = mflux6/(state7.rho * V7)
            state7_pitot = pitot_condition(state7, V7)
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
        state1 = make_gas_from_name(gasName)
        state1.set_pT(p1, T1)
        state0 = total_condition(state1, V1)
        fout.write('Total condition:\n')
        state0.write_state(fout)
    elif task in ['pitot', 'PITOT', 'Pitot']:
        fout.write('Input parameters:\n')
        fout.write('    Gas is %s, p1: %g Pa, T1: %g K, V1: %g m/s\n' % (gasName,p1,T1,V1) )
        state1 = make_gas_from_name(gasName)
        state1.set_pT(p1, T1)
        state0 = pitot_condition(state1, V1)
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
