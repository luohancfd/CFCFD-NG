#!/usr/bin/env python
"""
hadas85.py -- Hadas' 8.5 expansion-tube condition.

Done as an example of using gas_flow functions.
PJ, 21-Mar-2012
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.gas_flow import normal_shock, finite_wave_dp, normal_shock_p2p1
from cfpylib.nm.zero_solvers import secant

def main():
    print "Titan gas"
    state1 = Gas({'N2':0.95, 'CH4':0.05}, inputUnits='moles', outputUnits='moles')
    state1.set_pT(2600.0, 300.0)
    print "state1:"
    state1.write_state(sys.stdout)
    #
    print "Air accelerator gas"
    state10 = Gas({'Air':1.0})
    state10.set_pT(10.0, 300.0)
    print "state10:"
    state10.write_state(sys.stdout)
    #
    print "Incident shock"
    state2 = state1.clone()
    V2,V2g = normal_shock(state1, 4100.0, state2)
    print "V2=", V2, "Vg=", V2g, "expected 3670.56"
    print "state2:"
    state2.write_state(sys.stdout)
    print "Checks:"
    print "p2/p1=", state2.p/state1.p, "expected 166.4"
    print "rho2/rho1=", state2.rho/state1.rho, "expected 9.5474"
    print "T2/T1=", state2.T/state1.T, "expected 14.9"
    #
    print "\nNow do unsteady expansion..."
    # For the unsteady expansion of the test gas, regulation of the amount
    # of expansion is determined by the shock-processed accelerator gas.
    # Across the contact surface between these gases, the pressure and velocity
    # have to match so we set up some trials of various pressures and check 
    # that velocities match.
    def error_in_velocity(p5p2, state2=state2, V2g=V2g, state10=state10):
        "Compute the velocity mismatch for a given pressure ratio across the expansion."
        # Across the expansion, we get a test-gas velocity, V5g.
        V5g, state5 = finite_wave_dp('cplus', V2g, state2, p5p2*state2.p)
        # Across the contact surface, p20 == p5
        p20 = p5p2 * state2.p
        print "current guess for p5 and p20=", p20
        V10, V20, V20g, state20 = normal_shock_p2p1(state10, p20/state10.p)
        return (V5g - V10)/V5g # V10 was V20g - lab speed of accelerator gas - we now make the assumption that this is the same as the shock speed
    p5p2 = secant(error_in_velocity, 0.01, 0.011, tol=1.0e-3)
    print "From secant solve: p5/p2=", p5p2
    # It would have been faster and the code closer to Hadas' spreadsheet if we had
    # stepped down in pressure until we found the point where the velocities matched.
    # The expansion along the u+a wave would have appeared in the code here.
    V5g, state5 = finite_wave_dp('cplus', V2g, state2, p5p2*state2.p)
    print "Expanded test gas, at end of acceleration tube:"
    print "V5g=", V5g
    print "state5:"
    state5.write_state(sys.stdout)
    V10, V20, V20g, state20 = normal_shock_p2p1(state10, state5.p/state10.p)
    print V10
    print "Done."
    return

if __name__ == '__main__':
    main()

