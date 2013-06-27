# affine.py
# Scale the x-position of the CUBRC heat-transfer data
# using x_mm = alpha * x/L + beta 
# to match a two key points of the computational result.
# 1. separation point x/L=0.338 x=60mm
# 2. second-last transducer x/L=1.610 x=155mm
# Note that this transformation is not necessary for the pressure data.
# PJ, 25-June-2013
from numpy import array, linalg
a = array([[0.338, 1.0],[1.610, 1.0]])
b = array([60.0,155.0])
e = linalg.solve(a,b)
alpha, beta = e
print "alpha=", alpha, "beta=", beta
