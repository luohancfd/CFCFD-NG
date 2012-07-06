#!/usr/bin/env python
"""
FZ's theiving of PJ's script for Hadas.

Done as an example of using gas_flow functions.
PJ, 21-Mar-2012

FZ, 5th July 2012

The naming convention used here follow's Trimpi's report.
NASA TR R-133 'A preliminary theoretical study of the
expansion tube, a new device for producing high-enthalpy
short-duration hypersonic gas flows', Robert L. Trimpi
In addition the nozzle exit condition has been labelled
as 6, and the conehead surface condition as 7.
This is an X2 specific script configured for FZ's 8.4km/s
condition.
The script uses some experimental data as variables that we
try to match.
At this stage it is still somewhat manual but hopefully will
be updated in time.
The nozzle_exp_ratio needs to be manually tweaked until the
conehead pressure from the calculation matches that from the
experiment.
FZ.

"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.gas_flow import normal_shock, finite_wave_dp, normal_shock_p2p1, finite_wave_dv, pitot_condition, total_condition, steady_flow_with_area_change, theta_cone, beta_cone
from cfpylib.nm.zero_solvers import secant
import math

def main():
    # These are the inputs!
    primary_shock = 4100.0 # m/s
    secondary_shock = 8400.0 # m/s
    nozzle_exp_ratio = 3.0
    conehead_pressure = 8300.0 # Pa
    # Set up the shock tube test gas
    print "Air test gas"
    state1 = Gas({'Air':1.0}, inputUnits='moles', outputUnits='moles', with_ions=True)
    state1.set_pT(3000.0, 300.0)
    print "state1:"
    state1.write_state(sys.stdout)
    # Set up the acceleration tube gas
    print "Air accelerator gas"
    state10 = Gas({'Air':1.0})
    state10.set_pT(10.0, 300.0)
    print "state10:"
    state10.write_state(sys.stdout)
    # Compute the primary shock here
    print "Incident shock"
    state2 = state1.clone()
    V2,V2g = normal_shock(state1, primary_shock, state2)
    print "state2:"
    state2.write_state(sys.stdout)
    # For the unsteady expansion of the test gas, regulation of the amount
    # of expansion is determined by the shock-processed accelerator gas.
    print "\nNow do unsteady expansion..."
    V5g, state5 = finite_wave_dv('cplus', V2g, state2, secondary_shock)
    # Write out expanded test gas conditions
    print "\nExpanded test gas, at end of acceleration tube:"
    print "V5g=", V5g
    print "state5:"
    state5.write_state(sys.stdout)
    # Now lets put in the nozzle with a steady flow area change
    V6g, state6 = steady_flow_with_area_change(state5, V5g, nozzle_exp_ratio)
    # Print out nozzle expanded flow properties
    print "\nNozzle exit gas conditions for an area ratio of 3.0:"
    print "V6g:", V6g
    print "state6:"
    state6.write_state(sys.stdout)
    # Try and work out my conehead pressure for this particular flow condition
    # First I need to work out the shock angle
    shock_angle = beta_cone(state6, V6g, math.radians(15))
    print "\nShock angle over cone:", math.degrees(shock_angle)
    # Reverse the process to get the flow state behind the shock and check the surface angle is correct
    delta_s, vel_cond, state7 = theta_cone(state6, V6g, shock_angle)
    print "Surface angle should be the same.....: 15deg = ", math.degrees(delta_s), "deg"
    print "\nConehead surface conditions:"
    state7.write_state(sys.stdout)
    # Need to check whether the pressure are the same
    print "\nComputed conehead pressure should be the same as the experimental pressure"
    print "Computed:", state7.p, "Experimental:",conehead_pressure
    print "If pressure are not the same 'tweak' the nozzle_exp_ratio until they match"
    # Grab the total conditions
    total6 = total_condition(state6, V6g)
    # Now print my total conditions
    print "\nTotal conditions of nozzle exit test gas:"
    print "total6:"
    total6.write_state(sys.stdout)
    # And finally.... let's output a 'flight equivalent' speed
    flight_eq = math.sqrt(2*total6.h)
    print "\nFlight equivalent speed:", flight_eq
    print "Done."
    return

if __name__ == '__main__':
    main()

