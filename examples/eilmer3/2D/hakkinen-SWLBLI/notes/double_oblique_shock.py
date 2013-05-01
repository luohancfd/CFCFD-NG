# double_oblique_shock.py
"""
Estimate pressure rise across a reflected oblique shock.
PJ, 01-May-2013
"""
print "Begin..."
from cfpylib.gasdyn.ideal_gas_flow import *
import math
M1 = 2.0
p1 = 1.0
g = 1.4
print "First shock: ",
delta1 = 3.09 * math.pi/180.0
beta1 = beta_obl(M1,delta1,g)
p2 = p2_p1_obl(M1,beta1,g)
M2 = M2_obl(M1,beta1,delta1,g)
print "beta1=", beta1, "p2=", p2, "M2=", M2

print "Reflected shock:",
delta2 = delta1
beta2 = beta_obl(M2,delta2,g)
p3 = p2 * p2_p1_obl(M2,beta2,g)
M3 = M2_obl(M2,beta2,delta2,g)
print "beta2=", beta2, "p3=", p3, "M3=", M3

print "Done."
