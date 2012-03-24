#!/usr/bin/env python
"""
classic_shock_tube.py

Moderately high-performance shock tube with helium driving air.
Done as an example of using cea2_gas_flow functions but can be 
compared the Eilmer3 sod shock tube example.

PJ, 22-Mar-2012
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.cea2_gas_flow import shock_real, finite_wave_dp, shock_real_p2p1
from cfpylib.nm.zero_solvers import secant

def main():
    print "Helium driver gas"
    state4 = Gas({'He':1.0})
    state4.set_pT(30.0e6, 3000.0)
    print "state1:"
    state4.write_state(sys.stdout)
    #
    print "Air driven gas"
    state1 = Gas({'Air':1.0})
    state1.set_pT(30.0e3, 300.0)
    print "state1:"
    state1.write_state(sys.stdout)
    #
    print "\nNow do the classic shock tube solution..."
    # For the unsteady expansion of the driver gas, regulation of the amount
    # of expansion is determined by the shock-processed test gas.
    # Across the contact surface between these gases, the pressure and velocity
    # have to match so we set up some trials of various pressures and check 
    # that velocities match.
    def error_in_velocity(p3p4, state4=state4, state1=state1):
        "Compute the velocity mismatch for a given pressure ratio across the expansion."
        # Across the expansion, we get a test-gas velocity, V3g.
        p3 = p3p4*state4.p
        V3g, state3 = finite_wave_dp('cplus', 0.0, state4, p3)
        # Across the contact surface.
        p2 = p3
        print "current guess for p3 and p2=", p2
        V1s, V2, V2g, state2 = shock_real_p2p1(state1, p2/state1.p)
        return (V3g - V2g)/V3g
    p3p4 = secant(error_in_velocity, 0.1, 0.11, tol=1.0e-3)
    print "From secant solve: p3/p4=", p3p4
    print "Expanded driver gas:"
    p3 = p3p4*state4.p
    V3g, state3 = finite_wave_dp('cplus', 0.0, state4, p3)
    print "V3g=", V3g
    print "state3:"
    state3.write_state(sys.stdout)
    print "Shock-processed test gas:"
    V1s, V2, V2g, state2 = shock_real_p2p1(state1, p3/state1.p)
    print "V1s=", V1s, "V2g=", V2g
    print "state2:"
    state2.write_state(sys.stdout)
    return

if __name__ == '__main__':
    main()
    print "Done."

