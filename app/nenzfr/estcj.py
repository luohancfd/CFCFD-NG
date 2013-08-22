#! /usr/bin/env python 
"""
estcj.py: Equilibrium Shock Tube Conditions, Junior

This program can be used to estimate flow conditions 
for shock-processed flows typical of high-performance 
shock-tunnels and expansion tubes.
The gas is assumed to remain in thermochemical equilibrium
and the flow processing is done in decoupled quasi-one-dimensional
wave processes such as shock waves and expansion fans.
For the reflected shock tunnel, this means that the initial,
quiescent test gas is first processed by the incident shock and 
subsequently by the reflected shock.  
The incident shock sets the inflow conditions for the reflected shock 
but there is no further interaction.

The program can do a number of calculations:

* flow in a reflected shock tube with or without a nozzle
* pitot pressure from free-stream flow condition
* stagnation (total) condition from free-stream condition
* code surface condition from free-stream condition
 
When run as an application, this program takes its input as
command line arguments, performs the requested calculations and outputs
the gast-state results.
To see what specific inputs are required, start the program as::

$ estcj.py --help

Which particular input parameters you need to supply depends on the
chosen task, however, a typical flow condition for the T4 shock tunnel
with the Mach 4 nozzle may be computed using::

$ estcj.py --task=stn --gas=air --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.37e6 --ar=27.0

The full output is a bit too much to include here, but you should see that
this condition has an enthalpy of 5.43 MJ/kg and the nozzle-exit condition
has a pressure of 93.6 kPa and a static temperature of 1284 degrees K,
with a flow speed of 2.95 km/s.


Getting the program set up
--------------------------
estcj.py is not a stand-alone file.
It comes as part of the cfcfd3 compressible-flow collection and
depends upon functions from the cfpylib library to do the specific 
calculations.
The easiest way to get started is to build and install from the 
nenzfr directory where this source file resides::

$ cd app/nenzfr/
$ make install

You may then call upon estcj.py so long as you have suitable
enviroment variables set, as per the installation instructions
for Eilmer3.


Some History
------------
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
from cfpylib.gasdyn.cea2_gas import Gas, make_gas_from_name
from cfpylib.gasdyn.gas_flow import *

# ----------------------------------------------------------------------------

VERSION_STRING = "07-Mar-2012"
DEBUG_ESTCJ  = 0  # if 1: some detailed data is output to help debugging
PRINT_STATUS = 1  # if 1: the start of each stage of the computation is noted.

# ----------------------------------------------------------------------------

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
    H1 = state1.e + state1.p/state1.rho
    result = {'state1':state1, 'H1':H1}
    #
    if PRINT_STATUS: print 'Start incident-shock calculation.'
    state2 = make_gas_from_name(gasName)
    (V2,Vg) = normal_shock(state1, Vs, state2)
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
    # entropy is set, then pressure is relaxed via an isentropic process
    state5s.set_pT(state5.p, state5.T);  
    if pe==None:
        state5s.p = state5.p
    else:
        state5s.p = pe;
    state5s.EOS(problemType='ps')
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
            x7 = secant(error_at_exit, 0.001*x6, 0.00005*x6, tol=1.0e-4,
                        limits=[1.0/state5s.p,1.0])
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
            x7 = secant(error_at_exit, 0.001*x6, 0.00005*x6, tol=2.0e-4,
                        limits=[1.0/state5s.p,1.0])
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
                  choices=['st', 'stn', 'stnp', 'ishock', 'total', 'pitot', 'cone'],
                  help=("particular calculation to make: "
                        "st = reflected shock tube; "
                        "stn = reflected shock tube with nozzle; "
                        "stnp = reflected shock tube with nozzle expanded to pitot; "
                        "ishock = incident shock only; "
                        "total = free-stream to total condition; "
                        "pitot = free-stream to Pitot condition; "
                        "cone = free-stream to Taylor-Maccoll cone flow"))
    op.add_option('--gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'air11species', 'air13species', 
                           'n2', 'co2', 'h2ne', 'ar'],
                  help=("name of gas model: "
                        "air; " "air5species; " "air11species; " "air13species; " 
                        "n2; " "co2; " "h2ne; " "Ar"))
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
    op.add_option('--sigma-deg', dest='cone_half_angle_deg', type='float', default=None,
                  help=("half-angle of the cone, in degrees"))
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
    cone_half_angle_deg = opt.cone_half_angle_deg
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
    if V1 is None and task in ['pitot', 'total', 'cone']:
        print "Need to supply a free-stream velocity."
        bad_input = True
    if cone_half_angle_deg is None and task in ['cone',]:
        print "Need to supply a cone half-angle (in degrees)."
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
    elif task in ['cone', 'CONE', 'Cone']:
        fout.write('Input parameters:\n')
        fout.write('    Gas is %s, p1: %g Pa, T1: %g K, V1: %g m/s, sigma: %g degrees\n' 
                   % (gasName,p1,T1,V1,cone_half_angle_deg) )
        state1 = make_gas_from_name(gasName)
        state1.set_pT(p1, T1)
        fout.write('Free-stream condition:\n')
        state1.write_state(fout)
        cone_half_angle_rad = cone_half_angle_deg*math.pi/180.0
        beta_rad = beta_cone(state1, V1, cone_half_angle_rad)
        theta_c, V_cone_surface, state2 = theta_cone(state1, V1, beta_rad)
        assert abs(theta_c - cone_half_angle_rad) < 0.001
        fout.write('Shock angle: %g (rad), %g (deg)\n' % (beta_rad, beta_rad*180.0/math.pi))
        fout.write('Cone-surface velocity: %g m/s\n' % (V_cone_surface,))
        fout.write('Cone-surface condition:\n')
        state2.write_state(fout)
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
