# -*- coding: utf-8 -*-
"""
Expansion tube mess around

Chris James (c.james4@uq.edu.au)

based on lots of stuff done by various people PJ, Rowan, Richard Morgan.

currently built to similate tube with no secondary driver, easier

umar condition, pure He driver, 3kPa shock tube, 10Pa acceleration tube

NOT FINISHED, BUT GETTING SOMEWHERE

October 2012
"""

#this basic format pillaged from code by PJ in classic_expansion_tube.py - CJ

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.gas_flow import normal_shock, finite_wave_dp, normal_shock_p2p1
from cfpylib.gasdyn.gas_flow import expand_from_stagnation
from cfpylib.nm.zero_solvers import secant
from cfpylib.gasdyn.ideal_gas_flow import p0_p

def main():
    
    #build some gas objects    
    
    print "pure He driver gas"
    state4 = Gas({'He':1.0}, inputUnits='moles', outputUnits='moles')
    state4.set_pT(2.79e+07, 2700.0)
    print "state1:"
    state4.write_state(sys.stdout)    
     
    print "Air test gas"
    state1 = Gas({'Air':1.0,})
    state1.set_pT(3000.0, 300.0)
    print "state1:"
    state1.write_state(sys.stdout)
    #
    print "Air accelerator gas"
    state5 = Gas({'Air':1.0})
    state5.set_pT(10.0, 300.0)
    print "state10:"
    state5.write_state(sys.stdout)
    
    print "Steady expansion of driver"
    #for 100%He condition mach number terminating steady expansion is 2.15
    #need to work out the corresponding pressure for this to put into
    #the function
    M4dash = 2.15 #mach number terminating steady expansion
    (state4dash, V4dash) = expand_from_stagnation(1.0/(p0_p(M4dash,state4.gam)),state4)
    print "state4dash:"
    state4dash.write_state(sys.stdout)
    print "V4dash = {0} m/s".format(V4dash)
    
    #I think we need to start guessing shock speeds to be able to keep going here
    
    print "Start working on the shock into the test gas"
    print "Guess US1 of 5080 m/s"
    state2 = state1.clone()
    V2,V2g = normal_shock(state1, 5080.0, state2)    
    print "state2:"
    state2.write_state(sys.stdout)
    print "Checks:"
    print "p2/p1=", state2.p/state1.p
    print "rho2/rho1=", state2.rho/state1.rho
    print "T2/T1=", state2.T/state1.T
    print "V2g = {0} m/s".format(V2g)
    
    print "Unsteady expansion from state 4dash to state 3"
    print "Well, giving it a go for now anyway"
    (V3, state3) = finite_wave_dp('cplus', V4dash, state4dash, state2.p, steps=100)
    print "state3:"
    state3.write_state(sys.stdout)    
    print "V3 = {0} m/s".format(V3)
    
    print "Start working on the shock into the accelerator gas"
    print "Guess US2 of 11100 m/s"
    state6 = state5.clone()
    V6,V6g = normal_shock(state5, 11100.0, state6)    
    print "state6:"
    state6.write_state(sys.stdout)
    print "Checks:"
    print "p2/p1=", state6.p/state5.p
    print "rho2/rho1=", state6.rho/state5.rho
    print "T2/T1=", state6.T/state5.T
    print "V6g = {0} m/s".format(V6g)
    
    print "Unsteady expansion from state 2 to state 7"
    print "Well, giving it a go for now anyway"
    (V7, state7) = finite_wave_dp('cplus', V2g, state2, state6.p, steps=100)
    print "state7:"
    state7.write_state(sys.stdout)    
    print "V7 = {0} m/s".format(V7)
    M7 = V7/(state7.gam*state7.R*state7.T)**(0.5)
    print "M7 = {0}".format(M7)  
    
    """
    
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
    return"""

if __name__ == '__main__':
    main()


