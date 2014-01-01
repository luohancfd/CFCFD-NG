#!/usr/bin/env python
"""
oblique_shock_example.py

Demonstration of using the library functions to compute flow conditions
across an oblique shock in equilibrium air. 
Data are chosen to match examples from Hunt and Souders NASA-SP-3093.

PJ, 01-Jan-2014
"""
from math import pi
import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.gas_flow import theta_oblique, beta_oblique

print "Example 1: Hunt and Souders Table VIII, sub-table (j)"
print "Given shock angle, compute post-shock conditions."
s1 = Gas({'Air':1.0})
s1.set_pT(52.671, 268.858)
print "Initial gas state:"
s1.write_state(sys.stdout)
beta = 45.0 * pi/180 # shock angle
V1 = 7.9248e3
theta, V2, s2 = theta_oblique(s1, V1, beta)
print("Following oblique shock, beta=%g degrees, theta=%g degrees (Hunt&Souders 45 40.638)" 
      % (beta*180/pi, theta*180/pi))
s2.write_state(sys.stdout)
print "Across shock:"
print "p2/p1=%g, T2/T1=%g (Hunt&Souders: 376.84 21.206)" % (s2.p/s1.p, s2.T/s1.T)

print "\nExample 2: Hunt and Souders Table VI, sub-table (a)"
print "Given deflection angle, compute shock angle and then post-shock conditions."
s1.set_pT(3542.7, 219.428)
print "Initial gas state:"
s1.write_state(sys.stdout)
theta = 33.671 * pi/180 # shock angle
V1 = 1.8288e3
beta = beta_oblique(s1, V1, theta)
print("Following oblique shock, beta=%g degrees, theta=%g degrees (Hunt&Souders 45 33.671)" 
      % (beta*180/pi, theta*180/pi))
theta2, V2, s2 = theta_oblique(s1, V1, beta)
s2.write_state(sys.stdout)
print "Across shock:"
print "p2/p1=%g, T2/T1=%g (Hunt&Souders: 22.23 4.4472)" % (s2.p/s1.p, s2.T/s1.T)

print "Done."
