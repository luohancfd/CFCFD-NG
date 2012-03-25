#!/usr/bin/env python
"""
classic_shock_tube.py

Moderately high-performance shock tube with helium driving air.
Done as an example of using gas_flow functions but can be 
compared the Eilmer3 sod shock tube example.

PJ, 22-Mar-2012
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.gas_flow import normal_shock, finite_wave_dp, normal_shock_p2p1
from cfpylib.nm.zero_solvers import secant

def main():
    print "Helium driver gas"
    state4 = Gas({'He':1.0})
    state4.set_pT(30.0e6, 3000.0)
    print "state4:"
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
        V1s, V2, V2g, state2 = normal_shock_p2p1(state1, p2/state1.p)
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
    V1s, V2, V2g, state2 = normal_shock_p2p1(state1, p3/state1.p)
    print "V1s=", V1s, "V2g=", V2g
    print "state2:"
    state2.write_state(sys.stdout)
    assert abs(V2g - V3g)/V3g < 1.0e-3
    #
    # Make a record for plotting against the Eilmer3 simulation data.
    # We reconstruct the expected data along a tube 0.0 <= x <= 1.0
    # at t=100us, where the diaphragm is at x=0.5.
    x_centre = 0.5 # metres
    t = 100.0e-6 # seconds
    fp = open('exact.data', 'w')
    fp.write('# 1:x(m)  2:rho(kg/m**3) 3:p(Pa) 4:T(K) 5:V(m/s)\n')
    print 'Left end'
    x = 0.0
    fp.write('%g %g %g %g %g\n' % (x, state4.rho, state4.p, state4.T, 0.0))
    print 'Upstream head of the unsteady expansion.'
    x = x_centre - state4.a * t
    fp.write('%g %g %g %g %g\n' % (x, state4.rho, state4.p, state4.T, 0.0))
    print 'The unsteady expansion in n steps.'
    n = 100
    dp = (state3.p - state4.p) / n
    state = state4.clone()
    V = 0.0
    p = state4.p
    for i in range(n):
        rhoa = state.rho * state.a
        dV = -dp / rhoa
        V += dV
        p += dp
        state.set_ps(p, state4.s)
        x = x_centre + t * (V - state.a)
        fp.write('%g %g %g %g %g\n' % (x, state.rho, state.p, state.T, V))
    print 'Downstream tail of expansion.'
    x = x_centre + t * (V3g - state3.a)
    fp.write('%g %g %g %g %g\n' % (x, state3.rho, state3.p, state3.T, V3g))
    print 'Contact surface.'
    x = x_centre + t * V3g
    fp.write('%g %g %g %g %g\n' % (x, state3.rho, state3.p, state3.T, V3g))
    x = x_centre + t * V2g  # should not have moved
    fp.write('%g %g %g %g %g\n' % (x, state2.rho, state2.p, state2.T, V2g))
    print 'Shock front'
    x = x_centre + t * V1s  # should not have moved
    fp.write('%g %g %g %g %g\n' % (x, state2.rho, state2.p, state2.T, V2g))
    fp.write('%g %g %g %g %g\n' % (x, state1.rho, state1.p, state1.T, 0.0))
    print 'Right end'
    x = 1.0
    fp.write('%g %g %g %g %g\n' % (x, state1.rho, state1.p, state1.T, 0.0))
    fp.close()
    return

if __name__ == '__main__':
    main()
    print "Done."

