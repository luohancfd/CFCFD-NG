from cfpylib.gasdyn.ideal_gas_flow import beta_cone
from math import degrees, radians

beta = beta_cone(V1=1000.0, p1=95.84e3, T1=1103.0, theta=radians(20.0))
print "beta=", degrees(beta), "degrees"

